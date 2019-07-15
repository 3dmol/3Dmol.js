"""
A Flask-SocketIo based server for the active learning environment project
"""
from flask import Flask, request  # pylint: disable=import-error
from flask_socketio import SocketIO, send, emit, join_room, leave_room, close_room  # pylint: disable=import-error

# for SOCKETIO
import eventlet  # pylint: disable=import-error

eventlet.monkey_patch()

APP = Flask(__name__)
APP.config['SECRET_KEY'] = 'secret!'
SOCKETIO = SocketIO(APP, async_mode='eventlet')

SESSIONS = {}  # indexed by session name


@SOCKETIO.on('message')
def handle_message(msg):
    ''' Message Handler Function '''
    print('Message: ' + msg)
    send(msg, broadcast=True)


@SOCKETIO.on('check session name event')
def handle_check_name(json):
    ''' A function that checks whether the session name exists '''
    name = json['name']
    print("checking ", name, name in SESSIONS)
    if name in SESSIONS:
        emit('check session name response', "exists")
    else:
        emit('check session name response', "free")


@SOCKETIO.on('create session event')
def handle_create_session(json):
    ''' A function that creates a session '''
    name = json['name']
    if name in SESSIONS:
        print("Invalid name", name)
        emit('create session response', False)
    else:
        emit('create session response', True)
        print("Session Created", name)
        SESSIONS[name] = {'cnt': 0,
                          'sid': request.sid,  # id of the controller
                          'state': json['state'],
                          'view': json['view'],
                          'query_updates': {},  # each clients individual selections
                          'query_active': False
                          }
        join_room(name)
        emit('session count', 0, room=request.sid)


@SOCKETIO.on('join session event')
def handle_join_session(json):
    '''' Helps the user to join a session '''
    name = json['name']
    print("Join", name)
    if name not in SESSIONS:
        emit('join session response', "0")
    else:
        emit('join session response', "1")
        join_room(name)
        try:
            SESSIONS[name]['cnt'] += 1
        except AttributeError:
            emit('error: restart connection')
        emit('session count', SESSIONS[name]['cnt'], room=SESSIONS[name]['sid'])
        # send session state to connecting client
        if 'view' in SESSIONS[name]:
            emit('viewer view change response', SESSIONS[name]['view'], room=request.sid)
        if 'state' in SESSIONS[name]:
            emit('viewer state change response', SESSIONS[name]['state'], room=request.sid)
        print("Joined Session: " + str(json))

        if SESSIONS[name]['query_active']:
            emit('query start response', SESSIONS[name]['query_active'], room=request.sid)


@SOCKETIO.on('delete session event')
def handle_delete_session(json):
    ''' Helps the initiator to delete the session '''
    name = json['name']
    if name not in SESSIONS:
        return
    if request.sid != SESSIONS[name]['sid']:
        return  # only the controller can send updates
    emit('leave session response', room=name)
    close_room(name)
    print("Session Deleted:" + str(json))
    del SESSIONS[name]
    emit('delete session response')


@SOCKETIO.on('leave session event')
def handle_leave_session(json):
    ''' Helps the joined person(non initiator) to leave the session '''
    name = json['name']
    leave_room(name)
    if name not in SESSIONS:
        return
    SESSIONS[name]['cnt'] -= 1

    updates = SESSIONS[name]['query_updates']
    if request.sid in updates:
        del updates[request.sid]
        emit('query update response', len(updates), room=name)

    emit('session count', SESSIONS[name]['cnt'], room=SESSIONS[name]['sid'])
    print("Session Left: " + str(json))
    emit('leave session response')


@SOCKETIO.on('viewer view change event')
def handle_viewer_change(json):
    ''' Handles View Change '''
    name = json['name']
    viewstate = json['view']
    if name not in SESSIONS:
        return
    if request.sid != SESSIONS[name]['sid']:
        print("Incorrect sid for viewer change")
        return  # only the controller can send updates
    SESSIONS[name]['view'] = viewstate
    emit('viewer view change response', viewstate, room=name)


@SOCKETIO.on('viewer state change event')
def handle_state_change(json):
    ''' Handles Model State Change '''
    name = json['name']
    state = json['state']
    print("State change", name)
    if name not in SESSIONS:
        return
    if request.sid != SESSIONS[name]['sid']:
        print("Incorrect sid for state change")
        return  # only the controller can send updates
    SESSIONS[name]['state'] = state
    emit('viewer state change response', state, room=name)


@SOCKETIO.on('query start')
def handle_query_start(json):
    ''' Function to start a query '''
    name = json['name']
    print("query start", name)
    if name not in SESSIONS:
        return
    if request.sid != SESSIONS[name]['sid']:
        return  # only the controller can send updates
    SESSIONS[name]['query_updates'] = {}
    SESSIONS[name]['query_active'] = json['kind']

    emit('query start response', json['kind'], room=name)


@SOCKETIO.on('query end')
def handle_query_end(json):
    ''' Function to end a query '''
    name = json['name']
    if name not in SESSIONS:
        return
    if request.sid != SESSIONS[name]['sid']:
        return  # only the controller can send updates
    SESSIONS[name]['query_updates'] = {}
    SESSIONS[name]['query_active'] = False
    emit('query end response', '', room=name)


@SOCKETIO.on('query update')
def handle_query_update(json):
    ''' A single update from a client, respond with number unique responses '''
    name = json['name']
    sel = json['selected']
    if name not in SESSIONS:
        return
    # process update
    updates = SESSIONS[name]['query_updates']
    updates[request.sid] = sel

    emit('query update response', len(updates), room=SESSIONS[name]['sid'])


@SOCKETIO.on('query fetch')
def handle_query_fetch(json):
    ''' Return the merged results '''
    name = json['name']
    if name not in SESSIONS:
        return
    # process update
    updates = SESSIONS[name]['query_updates']
    cnts = {}
    for positions in updates.values():
        for pos in positions:
            pos = tuple(pos)
            if pos in cnts:
                cnts[pos] += 1
            else:
                cnts[pos] = 1
    cnts = list(cnts.items())
    emit('query fetch response', cnts, room=request.sid)


if __name__ == '__main__':
    SOCKETIO.run(APP)
