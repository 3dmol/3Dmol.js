# Active learning session server

This Flask/Socket.IO service synchronizes classroom viewers and collects atom
selections. The reusable browser client is `LearningSession` in the optional
`3Dmol.ui` bundle; it supports the standalone viewer and embedded classrooms.
`server.py` keeps sessions in memory, so restarting the server discards them.
Ordinary embedded viewers do not need this service.

For setup, viewer URLs, and proxy configuration, use the maintained
[active learning tutorial](../tutorials/learning_environment.md). Build the main
repository first, then run `python3 server.py` from this directory; `-p` selects
the port (default 5000). It binds to `127.0.0.1` by default; use `--host 0.0.0.0`
for direct connections from other machines. The `static` symlink exposes the
repository root, and `/` redirects to `/static/viewer.html`.

The standalone viewer connects to the page's origin for both `localhost` and `127.0.0.1`.
When Apache serves the page, proxy `/socket.io/` to this server, or set
`window.viewerSessionServer` explicitly for a direct connection. For local testing
without Apache, run `python3 server.py -p 5068` and open
`http://127.0.0.1:5068/static/viewer.html`. See the tutorial for proxy configuration
and port-conflict troubleshooting.

Embedded classrooms opt in with `data-learning` or the JavaScript UI
`learning` configuration. They default to `https://3dmol.org`; set
`data-learning-server` or `LearningConfig.server` to use this deployment.
Hidden instructor activation requires an explicit session identifier shared by
the intended viewers. See the tutorial's [embedded-viewer examples](../tutorials/learning_environment.md#embedded-viewers)
for the data attributes, hotkey, and JavaScript API.

## Dependencies

Use Python 3.10 or newer in a virtual environment. From this directory, install
the dependencies from the maintained requirements file:

```bash
python3 -m pip install --upgrade -r requirements.txt
```

The server uses Flask-SocketIO's native `threading` mode, `simple-websocket` for
WebSocket transport, and one Gunicorn `gthread` worker. This is a
[supported Flask-SocketIO deployment](https://flask-socketio.readthedocs.io/en/latest/deployment.html#gunicorn-web-server)
with both polling and WebSocket support. Eventlet and monkey patching are no
longer used. The explicit threading configuration also applies if eventlet remains
installed in an older environment; it can be uninstalled from the server's virtual
environment after upgrading.

`python3 server.py` starts Gunicorn on Linux or macOS, with 100 worker threads by
default. Use `--threads` to adjust capacity for your class; persistent WebSocket
connections occupy worker threads. Keep one worker because sessions live in
process memory. For deployments that invoke Gunicorn directly, the equivalent
command from this directory is:

```bash
gunicorn --worker-class gthread --workers 1 --threads 100 --bind 127.0.0.1:5000 server:app
```

Session event handling and expiry share a lock so concurrent connections cannot
overwrite revisions, create competing controllers, or delete a renewed session.

## Synchronization and cleanup

The edit menu supports volumetric maps and distance measurements alongside model
appearance, labels, and molecular surfaces. Shapes supplied in shared URLs are
also restored and synchronized. Students receive materialized scene
geometry, and their edit and trajectory controls are disabled while connected.
Pending local loads and surface edits are canceled when a student joins.

Scene updates are debounced and deduplicated. Compact snapshots separate model
appearance from coordinates and stored trajectory frames; subsequent messages
replace only changed top-level sections. The server owns revision numbers and
retains a complete snapshot for joins and reconnects. A revision mismatch requests
resynchronization. Camera messages do not resend the scene.

Deploy the server and rebuilt viewer assets together, and reload existing browser
tabs after upgrading so both ends use the revisioned synchronization protocol.
Hidden embedded viewers additionally use session discovery: waiting clients
watch a session identifier, and the server announces when an instructor creates
it. Watching does not create a classroom, count as a joined student, or refresh
the classroom's inactivity timer. An older server without discovery support
cannot automatically connect waiting hidden-mode viewers.

Closing the instructor's tab or window ends the session when the server detects
the disconnect. Connected students leave the session, their previous control
settings are restored, and the session name becomes available again. This cleanup
also covers lost instructor connections, even if the browser cannot send its
final deletion request. Students can disconnect and reconnect while the instructor
remains connected.

Sessions also expire after 48 hours without valid activity. The server checks expiry
on access and sweeps once per minute, removes room membership, and notifies
connected viewers. Name-availability checks and rejected updates do not keep a
session alive. Session state is process-local and is lost on restart; run one
session-server process unless shared storage and routing are added.
