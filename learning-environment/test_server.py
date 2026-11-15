"""Session protocol tests using real Flask-Socket.IO clients, without a network port.

Run with the server dependencies installed:
    python3 -m unittest discover -s learning-environment -p test_server.py
"""
import copy
from concurrent.futures import ThreadPoolExecutor
import importlib.util
from pathlib import Path
import sys
import threading
import time
import unittest
from unittest.mock import patch

spec = importlib.util.spec_from_file_location('learning_session_server', Path(__file__).with_name('server.py'))
server = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = server
spec.loader.exec_module(server)

VIEW = [0, 0, 0, -30, 0, 0, 0, 1]
SCENE = {'models': [{'atoms': [{'index': 0, 'elem': 'C', 'x': 0, 'y': 0, 'z': 0}], 'frames': []}],
         'shapes': [], 'surfaces': [], 'labels': [], 'frame': 0}


def concurrently(*actions):
    """Start requests together and propagate exceptions from every worker."""
    ready = threading.Barrier(len(actions))

    def run(action):
        ready.wait(timeout=5)
        action()

    with ThreadPoolExecutor(max_workers=len(actions)) as workers:
        futures = [workers.submit(run, action) for action in actions]
        for future in futures:
            future.result(timeout=5)


class SessionServerTests(unittest.TestCase):
    def setUp(self):
        self.now = 1000.0
        self.clients = []
        server.sessions.clear()
        server.session_watchers.clear()
        server._cleanup_started = False
        self.clock = patch.object(server, '_clock', lambda: self.now)
        self.clock.start()
        self.background = patch.object(server.socketio, 'start_background_task')
        self.start_background = self.background.start()

    def tearDown(self):
        for client in self.clients:
            if client.is_connected():
                client.disconnect()
        server.sessions.clear()
        server.session_watchers.clear()
        self.clock.stop()
        self.background.stop()

    def client(self):
        client = server.socketio.test_client(server.app)
        self.clients.append(client)
        return client

    def events(self, client, event=None):
        messages = client.get_received()
        if event is None:
            return messages
        return [message['args'][0] if message['args'] else None
                for message in messages if message['name'] == event]

    def create(self, name='Classroom'):
        teacher = self.client()
        teacher.emit('create session event', {'name': name, 'state': copy.deepcopy(SCENE), 'view': VIEW})
        secret = self.events(teacher, 'create session response')[0]
        self.assertGreater(secret, 0)
        return teacher, secret

    def join(self, name='Classroom'):
        student = self.client()
        student.emit('join session event', {'name': name})
        return student

    def test_home_serves_the_checkout_viewer_and_session_client(self):
        client = server.app.test_client()
        response = client.get('/')
        self.assertEqual(response.status_code, 302)
        self.assertEqual(response.headers['Location'], '/static/viewer.html')
        root = Path(__file__).resolve().parent.parent
        for file in ('viewer.html', 'learning-environment/sessions.js'):
            with self.subTest(file=file):
                response = client.get('/static/' + file)
                try:
                    self.assertEqual(response.status_code, 200)
                    self.assertEqual(response.data, (root / file).read_bytes())
                finally:
                    response.close()

    def test_concurrent_creation_grants_only_one_controller(self):
        contenders = [self.client(), self.client()]
        validate = server._valid_state

        def yield_after_name_check(value, partial=False):
            valid = validate(value, partial)
            # Release the GIL after the free-name check, before storing a session.
            # This exposes a check/write race without timing-based assertions.
            time.sleep(0.02)
            return valid

        with patch.object(server, '_valid_state', side_effect=yield_after_name_check):
            concurrently(*(lambda client=client: client.emit('create session event', {
                'name': 'Classroom', 'state': copy.deepcopy(SCENE), 'view': VIEW
            }) for client in contenders))

        responses = [self.events(client, 'create session response')[0] for client in contenders]
        self.assertEqual(sum(secret > 0 for secret in responses), 1)
        self.assertEqual(responses.count(0), 1)
        winner = contenders[next(index for index, secret in enumerate(responses) if secret > 0)]
        winner_sid = server.socketio.server.manager.sid_from_eio_sid(winner.eio_sid, '/')
        self.assertEqual(server.sessions['Classroom']['sid'], winner_sid)
        self.assertEqual(server.sessions['Classroom']['secret'], max(responses))
        self.assertEqual(set(server.socketio.server.manager.rooms['/']['Classroom']), {winner_sid})

    def test_waiting_embeds_discover_only_their_classroom_without_joining(self):
        waiting, other = self.client(), self.client()
        waiting.emit('watch session event', {'name': 'Classroom'})
        other.emit('watch session event', {'name': 'Another classroom'})
        self.assertEqual(server.sessions, {})
        self.assertEqual(self.events(waiting), [])
        teacher, _ = self.create()
        self.assertEqual(self.events(waiting, 'session available'), [{'name': 'Classroom'}])
        self.assertEqual(self.events(other), [])
        self.assertEqual(server.sessions['Classroom']['cnt'], 0)
        self.assertEqual(server.sessions['Classroom']['members'], set())
        teacher.emit('viewer state change event', {'name': 'Classroom', 'baseRevision': 0, 'patch': {'frame': 1}})
        self.assertEqual(self.events(waiting), [])
        waiting.emit('join session event', {'name': 'Classroom'})
        self.assertEqual(self.events(waiting, 'viewer state change response')[0]['state']['frame'], 1)
        self.assertEqual(server.sessions['Classroom']['cnt'], 1)

    def test_discovery_of_existing_sessions_does_not_extend_their_lifetime(self):
        teacher, _ = self.create()
        waiting = self.client()
        self.now += server.SESSION_TTL_SECONDS - 1
        waiting.emit('watch session event', {'name': 'Classroom'})
        self.assertEqual(self.events(waiting, 'session available'), [{'name': 'Classroom'}])
        self.assertEqual(server.sessions['Classroom']['last_activity'], 1000.0)
        self.now += 1
        waiting.emit('watch session event', {'name': 'Classroom'})
        self.assertEqual(self.events(waiting), [])
        self.assertNotIn('Classroom', server.sessions)
        self.assertEqual(self.events(teacher, 'leave session response'), [{'reason': 'expired'}])
        # Waiting clients remain subscribed when an instructor starts again.
        self.create()
        self.assertEqual(self.events(waiting, 'session available'), [{'name': 'Classroom'}])

    def test_watch_switch_unwatch_and_disconnect_release_discovery_subscriptions(self):
        waiting = self.client()
        sid = server.socketio.server.manager.sid_from_eio_sid(waiting.eio_sid, '/')
        waiting.emit('watch session event', {'name': 'Old'})
        waiting.emit('watch session event', {'name': 'New'})
        waiting.emit('unwatch session event', {'name': 'Old'})
        self.assertEqual(server.session_watchers, {sid: 'New'})
        self.create('Old')
        self.assertEqual(self.events(waiting), [])
        waiting.emit('unwatch session event', {'name': 'New'})
        self.assertEqual(server.session_watchers, {})
        self.create('New')
        self.assertEqual(self.events(waiting), [])
        waiting.emit('watch session event', {'name': 'Next'})
        waiting.disconnect()
        self.assertEqual(server.session_watchers, {})

    def test_concurrent_watch_and_create_deliver_one_availability_notification(self):
        teacher, waiting = self.client(), self.client()
        concurrently(
            lambda: waiting.emit('watch session event', {'name': 'Classroom'}),
            lambda: teacher.emit('create session event', {'name': 'Classroom', 'state': copy.deepcopy(SCENE), 'view': VIEW})
        )
        self.assertEqual(self.events(waiting, 'session available'), [{'name': 'Classroom'}])
        self.assertEqual(server.sessions['Classroom']['cnt'], 0)

    def test_concurrent_updates_cannot_accept_the_same_base_revision_twice(self):
        teacher, _ = self.create()
        student = self.join()
        self.events(teacher)
        self.events(student)
        revision_matches = server._revision_matches

        def yield_after_revision_check(message, session):
            matches = revision_matches(message, session)
            # Let another request run between validation and mutation unless the
            # handler owns the shared session lock throughout the transaction.
            time.sleep(0.02)
            return matches

        with patch.object(server, '_revision_matches', side_effect=yield_after_revision_check):
            concurrently(*(lambda frame=frame: teacher.emit('viewer state change event', {
                'name': 'Classroom', 'baseRevision': 0, 'patch': {'frame': frame}
            }) for frame in (1, 2)))

        packets = self.events(student, 'viewer state change response')
        self.assertEqual(len(packets), 1)
        self.assertEqual(packets[0]['baseRevision'], 0)
        self.assertEqual(packets[0]['revision'], 1)
        self.assertIn(packets[0]['patch']['frame'], (1, 2))
        self.assertEqual(server.sessions['Classroom']['revision'], 1)
        self.assertEqual(server.sessions['Classroom']['state']['frame'], packets[0]['patch']['frame'])
        replies = [(message['name'], message['args'][0]) for message in self.events(teacher)]
        self.assertEqual(replies, [('viewer state revision', {'revision': 1}),
                                  ('viewer state resync', {'revision': 1})])

    def test_expiration_is_inclusive_and_closes_the_room(self):
        teacher, _ = self.create()
        student = self.join()
        self.events(teacher)
        self.events(student)
        self.now += server.SESSION_TTL_SECONDS - 1
        self.assertEqual(server.expire_sessions(), [])
        self.now += 1
        self.assertEqual(server.expire_sessions(), ['Classroom'])
        self.assertEqual(server.sessions, {})
        self.assertEqual(self.events(teacher, 'leave session response'), [{'reason': 'expired'}])
        self.assertEqual(self.events(student, 'leave session response'), [{'reason': 'expired'}])
        self.assertNotIn('Classroom', server.socketio.server.manager.rooms.get('/', {}))
        self.assertEqual(server.expire_sessions(), [])

    def test_periodic_cleanup_starts_once_and_expires_without_requests(self):
        self.create()
        self.client()
        self.start_background.assert_called_once_with(server._cleanup_loop)
        self.now += server.SESSION_TTL_SECONDS
        with patch.object(server.socketio, 'sleep', side_effect=[None, StopIteration]) as sleep:
            with self.assertRaises(StopIteration):
                server._cleanup_loop()
        self.assertEqual(sleep.call_args_list[0].args, (server.SESSION_SWEEP_SECONDS,))
        self.assertNotIn('Classroom', server.sessions)

    def test_activity_extends_lifetime_but_name_probes_do_not(self):
        teacher, _ = self.create()
        self.now += server.SESSION_TTL_SECONDS - 1
        # A valid no-op is activity, even though it sends no scene/camera data.
        teacher.emit('viewer view change event', {'name': 'Classroom', 'view': VIEW})
        renewed = self.now
        self.assertEqual(server.sessions['Classroom']['last_activity'], renewed)
        self.now += server.SESSION_TTL_SECONDS - 1
        teacher.emit('check session name event', {'name': 'Classroom'})
        self.assertEqual(server.sessions['Classroom']['last_activity'], renewed)
        self.now += 1
        teacher.emit('check session name event', {'name': 'Classroom'})
        self.assertEqual(self.events(teacher, 'check session name response'), ['exists', 'free'])
        self.assertNotIn('Classroom', server.sessions)

    def test_each_valid_session_operation_records_activity(self):
        teacher, _ = self.create()
        student = self.join()
        actions = [
            (teacher, 'viewer view change event', {'view': [1, 2, 3, -30, 0, 0, 0, 1]}),
            (teacher, 'viewer state change event', {'patch': {'frame': 1}, 'baseRevision': 0}),
            (student, 'join session event', {}),
            (teacher, 'query start', {'kind': 'atoms'}),
            (student, 'query update', {'selected': [[1, 2, 3]]}),
            (teacher, 'query fetch', {}),
            (student, 'request state', {}),
            (teacher, 'query end', {}),
            (student, 'leave session event', {})
        ]
        for client, event, extra in actions:
            with self.subTest(event=event):
                self.now += 1
                client.emit(event, {'name': 'Classroom', **extra})
                self.assertEqual(server.sessions['Classroom']['last_activity'], self.now)

    def test_requests_cannot_revive_an_expired_session(self):
        teacher, secret = self.create()
        student = self.join()
        sid = self.events(student, 'join session response')[0]
        self.now += server.SESSION_TTL_SECONDS
        newcomer = self.join()
        self.assertEqual(self.events(newcomer, 'join session response'), ['0'])
        teacher.emit('viewer state change event', {'name': 'Classroom', 'state': SCENE})
        student.emit('reconnection', {'name': 'Classroom', 'secret': None, 'sid': sid})
        replacement = self.client()
        replacement.emit('reconnection', {'name': 'Classroom', 'secret': secret, 'sid': None})
        self.assertEqual(self.events(replacement, 'leave session response'), [{'reason': 'missing'}])
        self.assertNotIn('Classroom', server.sessions)
        # Expiration releases the name for a genuinely new session.
        self.create()
        self.assertEqual(server.sessions['Classroom']['revision'], 0)

    def test_student_membership_counts_disconnect_and_repeated_reconnection(self):
        teacher, _ = self.create()
        student = self.join()
        sid = self.events(student, 'join session response')[0]
        student.emit('join session event', {'name': 'Classroom'})
        self.assertEqual(server.sessions['Classroom']['cnt'], 1)
        teacher.emit('query start', {'name': 'Classroom', 'kind': 'atoms'})
        student.emit('query update', {'name': 'Classroom', 'selected': [[1, 2, 3]]})
        for _ in range(2):
            self.now += 1
            student.disconnect()
            self.assertEqual(server.sessions['Classroom']['cnt'], 0)
            student = self.client()
            student.emit('reconnection', {'name': 'Classroom', 'secret': None, 'sid': sid})
            response = self.events(student, 'reconnection response')[0]
            sid = response['sid']
            self.assertFalse(response['controller'])
            self.assertEqual(server.sessions['Classroom']['cnt'], 1)
            self.assertEqual(server.sessions['Classroom']['last_activity'], self.now)
            self.assertEqual(server.sessions['Classroom']['query_updates'], {sid: [[1, 2, 3]]})
        student.emit('leave session event', {'name': 'Classroom'})
        student.emit('leave session event', {'name': 'Classroom'})
        self.assertEqual(server.sessions['Classroom']['cnt'], 0)
        self.assertEqual(server.sessions['Classroom']['query_updates'], {})

    def test_browser_windows_sharing_cookies_keep_independent_roles(self):
        browser = server.app.test_client()
        with browser.session_transaction() as cookie_session:
            cookie_session['browser'] = 'shared-browser'
        self.assertIsNotNone(browser.get_cookie(server.app.config['SESSION_COOKIE_NAME']))
        teacher, student = [server.socketio.test_client(server.app, flask_test_client=browser)
                            for _ in range(2)]
        self.clients.extend([teacher, student])
        sid = lambda client: server.socketio.server.manager.sid_from_eio_sid(client.eio_sid, '/')
        teacher_sid, student_sid = sid(teacher), sid(student)
        self.assertNotEqual(teacher_sid, student_sid)

        teacher.emit('create session event', {'name': 'Classroom', 'state': copy.deepcopy(SCENE), 'view': VIEW})
        self.assertGreater(self.events(teacher, 'create session response')[0], 0)
        student.emit('join session event', {'name': 'Classroom'})
        self.assertEqual(self.events(student, 'join session response'), [student_sid])
        self.assertEqual([message['name'] for message in self.events(teacher)], ['session count'])
        self.assertEqual(server.sessions['Classroom']['sid'], teacher_sid)
        self.assertTrue(server.sessions['Classroom']['controller_connected'])
        self.assertEqual(server.sessions['Classroom']['members'], {student_sid})

        student.emit('viewer state change event', {
            'name': 'Classroom', 'baseRevision': 0, 'patch': {'frame': 99}
        })
        self.assertEqual(server.sessions['Classroom']['state'], SCENE)
        self.assertEqual(server.sessions['Classroom']['revision'], 0)
        self.assertEqual(self.events(teacher), [])
        teacher.emit('viewer state change event', {
            'name': 'Classroom', 'baseRevision': 0, 'patch': {'frame': 1}
        })
        self.assertEqual(self.events(student, 'viewer state change response'), [{
            'baseRevision': 0, 'revision': 1, 'patch': {'frame': 1}
        }])
        self.assertEqual(self.events(teacher, 'viewer state revision'), [{'revision': 1}])

    def test_delta_updates_keep_full_state_for_late_joiners_and_resync(self):
        teacher, _ = self.create()
        student = self.join()
        self.events(student)
        self.events(teacher)
        surfaces = [{'type': 1, 'style': {'opacity': 0.5}}]
        teacher.emit('viewer state change event', {
            'name': 'Classroom', 'patch': {'surfaces': surfaces}, 'baseRevision': 0, 'revision': 1
        })
        self.assertEqual(self.events(student, 'viewer state change response'), [{
            'patch': {'surfaces': surfaces}, 'baseRevision': 0, 'revision': 1
        }])
        self.assertEqual(self.events(teacher, 'viewer state revision'), [{'revision': 1}])
        self.assertEqual(server.sessions['Classroom']['state']['models'], SCENE['models'])
        late = self.join()
        snapshot = self.events(late, 'viewer state change response')[0]
        self.assertEqual(snapshot, {'state': {**SCENE, 'surfaces': surfaces}, 'revision': 1})
        student.emit('request state', {'name': 'Classroom'})
        self.assertEqual(self.events(student, 'viewer state change response'), [snapshot])

    def test_duplicate_updates_ack_without_broadcasting_or_ending_queries(self):
        teacher, _ = self.create()
        student = self.join()
        teacher.emit('query start', {'name': 'Classroom', 'kind': 'atoms'})
        self.events(student)
        self.events(teacher)
        for packet in [
            {'patch': {}, 'baseRevision': 0},
            {'patch': {'frame': 0}, 'baseRevision': 0, 'revision': 1},
            {'state': SCENE, 'baseRevision': 0},
            {'state': SCENE}
        ]:
            teacher.emit('viewer state change event', {'name': 'Classroom', **packet})
        teacher.emit('viewer view change event', {'name': 'Classroom', 'view': VIEW})
        self.assertEqual(self.events(student), [])
        self.assertEqual(self.events(teacher, 'viewer state revision'), [{'revision': 0}] * 4)
        self.assertEqual(server.sessions['Classroom']['query_active'], 'atoms')
        self.assertEqual(server.sessions['Classroom']['revision'], 0)

    def test_appearance_updates_preserve_queries_but_frame_changes_end_them(self):
        teacher, _ = self.create()
        student = self.join()
        teacher.emit('query start', {'name': 'Classroom', 'kind': 'atoms'})
        student.emit('query update', {'name': 'Classroom', 'selected': [[0, 0, 0]]})
        self.events(student)
        teacher.emit('viewer state change event', {
            'name': 'Classroom', 'baseRevision': 0,
            'patch': {'labels': [{'text': 'Instructor annotation'}],
                      'surfaces': [{'style': {'opacity': 0.5}}],
                      'modelStyles': [[{'sphere': {'color': 'blue'}}]]}
        })
        self.assertEqual(self.events(student, 'query end response'), [])
        self.assertEqual(server.sessions['Classroom']['query_active'], 'atoms')
        self.assertEqual(len(server.sessions['Classroom']['query_updates']), 1)
        teacher.emit('viewer state change event', {
            'name': 'Classroom', 'baseRevision': 1, 'patch': {'frame': 1}
        })
        self.assertEqual(self.events(student, 'query end response'), [''])
        self.assertFalse(server.sessions['Classroom']['query_active'])
        self.assertEqual(server.sessions['Classroom']['query_updates'], {})

    def test_stale_or_malformed_updates_do_not_overwrite_state_or_extend_lifetime(self):
        teacher, _ = self.create()
        teacher.emit('viewer state change event', {'name': 'Classroom', 'patch': {'frame': 1}, 'baseRevision': 0})
        self.events(teacher)
        recorded = server.sessions['Classroom']['last_activity']
        self.now += 10
        teacher.emit('viewer state change event', {'name': 'Classroom', 'patch': {'frame': 99}, 'baseRevision': 0})
        self.assertEqual(self.events(teacher, 'viewer state resync'), [{'revision': 1}])
        for packet in [None, [], {}, {'state': []}, {'state': {}}, {'patch': []},
                       {'patch': {'models': None}, 'baseRevision': 1},
                       {'patch': {'surfaces': None}, 'baseRevision': 1},
                       {'patch': {'shapes': {}}, 'baseRevision': 1},
                       {'patch': {'labels': 'invalid'}, 'baseRevision': 1},
                       {'patch': {'modelStyles': {}}, 'baseRevision': 1},
                       {'patch': {'view': []}, 'baseRevision': 1},
                       {'patch': {'frame': -1}, 'baseRevision': 1},
                       {'state': SCENE, 'patch': {}},
                       {'patch': {'frame': 99}, 'baseRevision': None},
                       {'patch': {'frame': 99}, 'baseRevision': True},
                       {'patch': {'frame': 99}, 'baseRevision': 1, 'revision': 8}]:
            teacher.emit('viewer state change event', {'name': 'Classroom', **packet} if isinstance(packet, dict) else packet)
        for value in [None, [], VIEW[:-1], [float('nan')] * 8]:
            teacher.emit('viewer view change event', {'name': 'Classroom', 'view': value})
        self.assertEqual(server.sessions['Classroom']['state']['frame'], 1)
        self.assertEqual(server.sessions['Classroom']['revision'], 1)
        self.assertEqual(server.sessions['Classroom']['last_activity'], recorded)

    def test_only_instructor_can_edit_or_fetch_query_results(self):
        teacher, _ = self.create()
        student = self.join()
        outsider = self.client()
        self.events(teacher)
        self.events(student)
        recorded = server.sessions['Classroom']['last_activity']
        self.now += 20
        for client in [student, outsider]:
            client.emit('viewer state change event', {'name': 'Classroom', 'state': {**SCENE, 'frame': 50}})
            client.emit('viewer view change event', {'name': 'Classroom', 'view': [5] * 8})
            client.emit('query start', {'name': 'Classroom', 'kind': 'atoms'})
            client.emit('query fetch', {'name': 'Classroom'})
            client.emit('delete session event', {'name': 'Classroom'})
        outsider.emit('request state', {'name': 'Classroom'})
        outsider.emit('query update', {'name': 'Classroom', 'selected': [[0, 0, 0]]})
        self.assertEqual(server.sessions['Classroom']['state'], SCENE)
        self.assertEqual(server.sessions['Classroom']['last_activity'], recorded)
        self.assertEqual(self.events(teacher), [])
        self.assertEqual(self.events(student), [])
        self.assertEqual(self.events(outsider), [])

    def test_controller_handoff_preserves_revision_and_validates_secret(self):
        teacher, secret = self.create()
        teacher.emit('viewer state change event', {'name': 'Classroom', 'state': {**SCENE, 'frame': 2}})
        bad = self.client()
        bad.emit('reconnection', {'name': 'Classroom', 'secret': secret + 1, 'sid': None})
        self.assertEqual(self.events(bad, 'leave session response'), [{'reason': 'invalid controller'}])
        replacement = self.client()
        replacement.emit('reconnection', {'name': 'Classroom', 'secret': secret, 'sid': None})
        response = self.events(replacement, 'reconnection response')[0]
        self.assertTrue(response['controller'])
        self.assertEqual(response['revision'], 1)
        teacher.disconnect()
        self.assertEqual(server.sessions['Classroom']['sid'], response['sid'])
        self.assertEqual(self.events(replacement), [])
        replacement.emit('viewer state change event', {'name': 'Classroom', 'patch': {'frame': 3}, 'baseRevision': 1})
        self.assertEqual(server.sessions['Classroom']['state']['frame'], 3)
        self.assertEqual(server.sessions['Classroom']['revision'], 2)

    def test_instructor_disconnect_removes_classroom_without_an_unload_event(self):
        teacher, secret = self.create()
        student = self.join()
        student_sid = self.events(student, 'join session response')[0]
        other_teacher, _ = self.create('Other classroom')
        other_session = server.sessions['Other classroom']

        # Closing a browser can drop its explicit delete-session message.
        teacher.disconnect()
        self.assertNotIn('Classroom', server.sessions)
        self.assertNotIn('Classroom', server.socketio.server.manager.rooms.get('/', {}))
        self.assertEqual(self.events(student, 'leave session response'),
                         [{'reason': 'instructor disconnected'}])
        self.assertTrue(student.is_connected())
        self.assertIs(server.sessions['Other classroom'], other_session)
        self.assertEqual(self.events(other_teacher), [])

        replacement = self.client()
        replacement.emit('reconnection', {'name': 'Classroom', 'secret': secret, 'sid': None})
        self.assertEqual(self.events(replacement, 'leave session response'), [{'reason': 'missing'}])
        student.emit('reconnection', {'name': 'Classroom', 'secret': None, 'sid': student_sid})
        self.assertEqual(self.events(student, 'leave session response'), [{'reason': 'missing'}])
        self.assertNotIn('Classroom', server.sessions)
        replacement.emit('check session name event', {'name': 'Classroom'})
        self.assertEqual(self.events(replacement, 'check session name response'), ['free'])
        self.create()
        self.assertEqual(server.sessions['Classroom']['revision'], 0)
        self.assertEqual(server.sessions['Classroom']['members'], set())
        self.assertIs(server.sessions['Other classroom'], other_session)

    def test_deleting_a_session_notifies_members_and_releases_state(self):
        teacher, _ = self.create()
        student = self.join()
        self.events(student)
        teacher.emit('delete session event', {'name': 'Classroom'})
        self.assertNotIn('Classroom', server.sessions)
        self.assertEqual(self.events(student, 'leave session response'), [{'reason': 'deleted'}])
        self.assertNotIn('Classroom', server.socketio.server.manager.rooms.get('/', {}))


if __name__ == '__main__':
    unittest.main()
