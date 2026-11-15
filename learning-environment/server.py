from flask import Flask, request, redirect, url_for
from flask_socketio import SocketIO, send, emit, join_room, leave_room
from functools import wraps
import argparse
import datetime
import math
import random
import threading
import time

app = Flask(__name__)
app.config['SECRET_KEY'] = 'secret!'
# Explicit selection prevents an old eventlet installation from taking priority.
# Handle each connection's events in receive order; other clients run in threads.
socketio = SocketIO(app, async_mode='threading', async_handlers=False,
                    cors_allowed_origins="*", max_http_buffer_size=1e8)

SESSION_TTL_SECONDS = 48 * 60 * 60
SESSION_SWEEP_SECONDS = 60
sessions = {}  # In-memory state, indexed by session name.
_cleanup_started = False
_cleanup_lock = threading.Lock()
_sessions_lock = threading.RLock()


def _serialized(function):
    """Keep session checks, mutations and emitted revisions atomic across threads."""
    @wraps(function)
    def locked(*args, **kwargs):
        with _sessions_lock:
            return function(*args, **kwargs)
    return locked


def _clock():
    # Elapsed inactivity must not change when the wall clock is corrected.
    return time.monotonic()


def _name(message):
    name = message.get('name') if isinstance(message, dict) else None
    return name if isinstance(name, str) and name.strip() else None


def _remove_session(name, reason):
    if sessions.pop(name, None) is not None:
        # The SocketIO instance methods also work in the background sweeper,
        # outside a Flask request context.
        socketio.emit('leave session response', {'reason': reason}, to=name, namespace='/')
        socketio.close_room(name, namespace='/')


@_serialized
def expire_sessions(now=None):
    """Remove sessions inactive for at least 48 hours; return the removed names."""
    now = _clock() if now is None else now
    expired = []
    for name, session in list(sessions.items()):
        if now - session['last_activity'] >= SESSION_TTL_SECONDS:
            _remove_session(name, 'expired')
            expired.append(name)
    return expired


def _cleanup_loop():
    while True:
        socketio.sleep(SESSION_SWEEP_SECONDS)
        expire_sessions()


def start_session_cleanup():
    """Start one sweeper per server process, including imported WSGI deployments."""
    global _cleanup_started
    with _cleanup_lock:
        if not _cleanup_started:
            _cleanup_started = True
            socketio.start_background_task(_cleanup_loop)


def _session(message):
    name = _name(message)
    session = sessions.get(name)
    if session is not None and _clock() - session['last_activity'] >= SESSION_TTL_SECONDS:
        _remove_session(name, 'expired')
        session = None
    return name, session


def _touch(session):
    session['last_activity'] = _clock()


def _controller(session):
    return session is not None and session['controller_connected'] and request.sid == session['sid']


def _member(session):
    return session is not None and (_controller(session) or request.sid in session['members'])


def _count(session):
    session['cnt'] = len(session['members'])
    if session['controller_connected']:
        emit('session count', session['cnt'], to=session['sid'])


def _send_state(session):
    emit('viewer view change response', session['view'], to=request.sid)
    emit('viewer state change response', {
        'state': session['state'], 'revision': session['revision']
    }, to=request.sid)


def _valid_state(value, partial=False):
    if not isinstance(value, dict) or any(not isinstance(key, str) for key in value):
        return False
    if not partial and not isinstance(value.get('models'), list):
        return False
    for key in ('models', 'modelStyles', 'shapes', 'surfaces', 'labels', 'modelTitles'):
        if key in value and not isinstance(value[key], list):
            return False
    for key in ('view', 'viewSettings', 'trajectory'):
        if key in value and not isinstance(value[key], dict):
            return False
    if 'frame' in value and (type(value['frame']) is not int or value['frame'] < 0):
        return False
    return True


def _valid_view(value):
    return isinstance(value, list) and len(value) == 8 and all(
        isinstance(item, (int, float)) and not isinstance(item, bool) and math.isfinite(item)
        for item in value
    )


def _revision_matches(message, session):
    base = message.get('baseRevision')
    # Legacy clients send complete snapshots without revisions.
    if 'baseRevision' not in message and 'state' in message and 'revision' not in message:
        return True
    if type(base) is not int or base < 0:
        return False
    revision = message.get('revision', base + 1)
    return type(revision) is int and revision == base + 1 and base == session['revision']


@socketio.on('connect')
@_serialized
def handleConnect():
    start_session_cleanup()


@socketio.on('message')
@_serialized
def handleMessage(message):
    print(datetime.datetime.now(), 'Message:', message)
    send(message, broadcast=True)


@socketio.on('check session name event')
@_serialized
def handleCheckName(message):
    name, session = _session(message)
    if name is not None:
        # Name probes do not extend the lifetime of an existing session.
        emit('check session name response', 'exists' if session is not None else 'free')


@socketio.on('create session event')
@_serialized
def handleCreateSession(message):
    name, existing = _session(message)
    if (name is None or existing is not None or not _valid_state(message.get('state'))
            or not _valid_view(message.get('view'))):
        emit('create session response', 0)
        return
    secret = random.randint(1, int(4e9))
    session = sessions[name] = {
        'cnt': 0, 'sid': request.sid, 'controller_connected': True,
        'members': set(), 'disconnected_members': set(),
        'state': message['state'], 'view': message['view'], 'revision': 0,
        'query_updates': {}, 'query_active': False, 'secret': secret,
        'last_activity': _clock()
    }
    join_room(name)
    emit('create session response', secret)
    emit('viewer state revision', {'revision': 0})
    _count(session)
    print('Session Created', name, request.sid)


@socketio.on('reconnection')
@_serialized
def handleReconnection(message):
    name, session = _session(message)
    if session is None:
        emit('leave session response', {'reason': 'missing'})
        return
    secret = message.get('secret')
    if secret is not None:
        if secret != session['secret']:
            emit('leave session response', {'reason': 'invalid controller'})
            return
        session['sid'] = request.sid
        session['controller_connected'] = True
        session['members'].discard(request.sid)
        controller = True
    else:
        oldsid = message.get('sid')
        if not isinstance(oldsid, str) or (oldsid not in session['members'] and oldsid not in session['disconnected_members']):
            emit('leave session response', {'reason': 'unknown member'})
            return
        session['members'].discard(oldsid)
        session['disconnected_members'].discard(oldsid)
        session['members'].add(request.sid)
        if oldsid in session['query_updates']:
            session['query_updates'][request.sid] = session['query_updates'].pop(oldsid)
        controller = False
    join_room(name)
    _touch(session)
    _count(session)
    emit('reconnection response', {'sid': request.sid, 'controller': controller,
                                   'revision': session['revision']})
    if not controller:
        _send_state(session)
        if session['query_active']:
            emit('query start response', session['query_active'])


@socketio.on('join session event')
@_serialized
def handleJoinSession(message):
    name, session = _session(message)
    if session is None:
        emit('join session response', '0')
        return
    if not _controller(session):
        session['members'].add(request.sid)
    session['disconnected_members'].discard(request.sid)
    _touch(session)
    join_room(name)
    emit('join session response', request.sid)
    _count(session)
    _send_state(session)
    if session['query_active']:
        emit('query start response', session['query_active'])


@socketio.on('delete session event')
@_serialized
def handleDeleteSession(message):
    name, session = _session(message)
    if _controller(session):
        _remove_session(name, 'deleted')
        emit('delete session response')


@socketio.on('leave session event')
@_serialized
def handleLeaveSession(message):
    name, session = _session(message)
    if session is None or request.sid not in session['members']:
        return
    session['members'].remove(request.sid)
    session['disconnected_members'].discard(request.sid)
    leave_room(name)
    _touch(session)
    if session['query_updates'].pop(request.sid, None) is not None:
        emit('query update response', len(session['query_updates']), to=session['sid'])
    _count(session)
    emit('leave session response', {'reason': 'left'})


@socketio.on('disconnect')
@_serialized
def handleDisconnect(reason=None):
    for name, session in list(sessions.items()):
        if request.sid == session['sid']:
            # Browser unload messages are unreliable. The classroom belongs to
            # its current instructor connection, including after a handoff.
            _remove_session(name, 'instructor disconnected')
            continue
        if request.sid in session['members']:
            session['members'].remove(request.sid)
            session['disconnected_members'].add(request.sid)
            _count(session)
    # Student disconnects preserve membership for reconnection without
    # extending the classroom's inactivity timeout.


@socketio.on('viewer view change event')
@_serialized
def handleViewerChange(message):
    name, session = _session(message)
    if not _controller(session) or not _valid_view(message.get('view')):
        return
    _touch(session)
    if message['view'] != session['view']:
        session['view'] = message['view']
        emit('viewer view change response', session['view'], to=name, skip_sid=request.sid)


@socketio.on('viewer state change event')
@_serialized
def handleStateChange(message):
    name, session = _session(message)
    if not _controller(session):
        return
    full = 'state' in message
    if full == ('patch' in message):
        return  # Require exactly one representation, never a malformed full reset.
    incoming = message['state'] if full else message['patch']
    if not _valid_state(incoming, partial=not full):
        return
    if not _revision_matches(message, session):
        emit('viewer state resync', {'revision': session['revision']})
        return
    _touch(session)
    if full:
        changed = incoming != session['state']
        updated = incoming
    else:
        incoming = {key: value for key, value in incoming.items()
                    if key not in session['state'] or value != session['state'][key]}
        changed = bool(incoming)
        updated = {**session['state'], **incoming}
    if changed:
        invalidates_query = (session['query_active'] or session['query_updates']) and (
            updated.get('models') != session['state'].get('models')
            or updated.get('frame', 0) != session['state'].get('frame', 0)
        )
        previous_revision = session['revision']
        session['revision'] += 1
        session['state'] = updated
        packet = {'state': updated, 'revision': session['revision']} if full else {
            'patch': incoming, 'baseRevision': previous_revision, 'revision': session['revision']
        }
        emit('viewer state change response', packet, to=name, skip_sid=request.sid)
        # Atom/frame changes invalidate positions; appearance edits do not.
        if invalidates_query:
            session['query_updates'] = {}
            session['query_active'] = False
            emit('query end response', '', to=name)
    # Acknowledge unchanged snapshots too, so a sender never advances past us.
    emit('viewer state revision', {'revision': session['revision']})


@socketio.on('request state')
@_serialized
def handleRequestState(message):
    _, session = _session(message)
    if _member(session):
        _touch(session)
        _send_state(session)


@socketio.on('query start')
@_serialized
def handleQueryStart(message):
    name, session = _session(message)
    kind = message.get('kind') if isinstance(message, dict) else None
    if not _controller(session) or not isinstance(kind, str) or not kind:
        return
    _touch(session)
    session['query_updates'] = {}
    session['query_active'] = kind
    emit('query start response', kind, to=name)


@socketio.on('query end')
@_serialized
def handleQueryEnd(message):
    name, session = _session(message)
    if not _controller(session):
        return
    _touch(session)
    if session['query_active'] or session['query_updates']:
        session['query_updates'] = {}
        session['query_active'] = False
        emit('query end response', '', to=name)


@socketio.on('query update')
@_serialized
def handleQueryUpdate(message):
    _, session = _session(message)
    selected = message.get('selected') if isinstance(message, dict) else None
    if (session is None or request.sid not in session['members'] or not session['query_active']
            or not isinstance(selected, list) or not all(
                isinstance(position, list) and len(position) == 3 and all(
                    isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value)
                    for value in position) for position in selected)):
        return
    _touch(session)
    if session['query_updates'].get(request.sid) != selected:
        session['query_updates'][request.sid] = selected
        emit('query update response', len(session['query_updates']), to=session['sid'])


@socketio.on('query fetch')
@_serialized
def handleQueryFetch(message):
    _, session = _session(message)
    if not _controller(session):
        return
    _touch(session)
    counts = {}
    for positions in session['query_updates'].values():
        for position in positions:
            position = tuple(position)
            counts[position] = counts.get(position, 0) + 1
    emit('query fetch response', list(counts.items()))


@app.route('/')
def home():
    return redirect(url_for('static', filename='viewer.html'))


def run_server(host, port, threads):
    """Serve HTTP, polling and WebSockets with one threaded Gunicorn worker."""
    from gunicorn.app.base import BaseApplication

    class SessionApplication(BaseApplication):
        def load_config(self):
            address = f'[{host}]' if ':' in host else host
            for key, value in {'bind': f'{address}:{port}', 'workers': 1,
                               'worker_class': 'gthread', 'threads': threads}.items():
                self.cfg.set(key, value)

        def load(self):
            # Start inside the worker, after Gunicorn forks; never in the parent.
            start_session_cleanup()
            return app

    SessionApplication().run()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run 3Dmol.js learning environment server.')
    parser.add_argument('-p', '--port', type=int, default=5000, help='Port to run on (default 5000)')
    parser.add_argument('--host', default='127.0.0.1', help='Bind address (default 127.0.0.1)')
    parser.add_argument('--threads', type=int, default=100, help='Worker threads (default 100)')
    args = parser.parse_args()
    if not 1 <= args.port <= 65535:
        parser.error('--port must be between 1 and 65535')
    if args.threads < 2:
        parser.error('--threads must be at least 2 to support concurrent Socket.IO requests')
    run_server(args.host, args.port, args.threads)
