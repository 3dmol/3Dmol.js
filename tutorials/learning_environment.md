<script src="../build/3Dmol-min.js"></script> 

The [hosted 3Dmol.js viewer](tutorial-url.html) and embedded viewers can serve as
an active learning environment for teaching properties of molecular structures.
The instructor creates a session using the **Classroom** button in the upper
right corner. Students who join see the instructor's molecular scene. They can
rotate and zoom independently, and instructor camera changes update their views.

The instructor can pose questions that the students then answer by clicking on atoms.  The accumulated and anonymous answers can then be displayed in the instructor's view.

## Embedded viewers

Classrooms are also available in embedded viewers through the optional
`3Dmol.ui` bundle. Include it after the core bundle, then opt in with
`data-learning`. No classroom connection is made by ordinary embeds.

For the same classroom button and session-name dialog as the standalone viewer:

```html
<meta charset="utf-8">
<script src="https://3dmol.org/build/3Dmol-min.js"></script>
<script src="https://3dmol.org/build/3Dmol.ui-min.js"></script>

<div id="classroom-viewer" class="viewer_3Dmoljs"
     style="width: 100%; height: 400px; position: relative;"
     data-pdb="1YCR" data-style="cartoon"
     data-ui="true" data-learning="true"></div>
```

Embedded classrooms use `https://3dmol.csb.pitt.edu` as their session server by default.
Set `data-learning-server="https://classroom.example.edu"` to use another server.
The server address includes the protocol and, when needed, a port; it is separate
from the molecular-data URL. All participants in a classroom must use the same
server and session identifier. A page can contain several independent classroom
viewers, each with its own identifier and controls.

For a standalone HTML file, put `<meta charset="utf-8">` near the beginning of
`<head>` and save the file as UTF-8. This also applies when opening the page
directly with `file://`, so symbols and other non-ASCII text display correctly.
Local pages still need an internet connection for the hosted scripts, molecule
downloads, and classroom server.
Use the direct server address `https://3dmol.csb.pitt.edu` when setting
`data-learning-server`; the `3dmol.org` alias redirects and can prevent a
WebSocket upgrade. Existing bundles can use the same attribute to avoid that
redirect.

### Hidden instructor activation

Use `data-learning="hidden"` when the classroom button should stay hidden until
an instructor starts the session. **An explicit `data-learning-session` is
required.** Choose an identifier unique to the class and viewer, and put the same
identifier on every participant's embed. Participants may be on different pages.

```html
<script src="https://3Dmol.csb.pitt.edu/build/3Dmol.js"></script>
<script src="https://3dmol.csb.pitt.edu/build/3Dmol.ui-min.js"></script>
<div id="bonding-viewer" class="viewer_3Dmoljs"
     style="width: 100%; height: 400px; position: relative;"
     data-pdb="1YCR" data-style="cartoon"
     data-learning="hidden"
     data-learning-session="biochem-bonding-class-a"></div>
```

This example uses the same core and UI scripts as the button example. Focus the
viewer by clicking it, then press **Ctrl+Alt+I** once to activate the instructor.
The classroom dialog immediately shows `Starting classroom...` while connecting.
The instructor gets the usual editing and classroom controls, including **Query
Atoms**, response counts, and ending the session. Use `data-learning-hotkey` to
choose another shortcut. Keyboard input in text fields does not activate the
shortcut.

If another instructor already holds the name when activation is requested, the
dialog reports that the name is in use. End the existing session or choose a
different session identifier before retrying. A failed activation keeps this
viewer out of student mode.

Other active hidden-mode embeds with the same server and identifier automatically
join as students when the instructor starts. Viewers opened after the instructor
starts also join. Student editing and trajectory controls are disabled while
following; students can still rotate, zoom, and answer atom queries. Waiting for
an instructor does not create a session or extend an existing session's lifetime.

Hidden mode controls how the interface appears; it does not authenticate an
instructor. Someone with the session identifier and shortcut can start a session
when that identifier is free.

### JavaScript UI API

The same feature is available through {@link UIConfig} and {@link LearningSession}
in `3Dmol.ui`. Enable it when creating the viewer's {@link StateManager}:

```javascript
const controls = new $3Dmol.StateManager(viewer, {
    edit: false,
    trajectory: false,
    learning: {
        mode: 'hidden',
        session: 'biochem-bonding-class-a',
        server: 'https://3dmol.csb.pitt.edu',
        hotkey: 'Ctrl+Alt+I'
    }
});
const classroom = controls.learning;
```

Use `learning: true` for the classroom button, or pass a {@link LearningConfig}
object. `mode` defaults to `'button'`, and `server` defaults to
`'https://3dmol.csb.pitt.edu'`. Hidden mode requires `session`. In button mode, `session`
can prefill the session-name field. `edit` and `trajectory` retain their ordinary
UI defaults; the example disables them until this viewer becomes the instructor.

`classroom.open()` opens the classroom dialog. Applications can also call
`classroom.start(name)` to create an instructor session or `classroom.join(name)`
to join an existing session as a student. The configured `session` can be used by
omitting `name`. `classroom.getRole()` returns `'idle'`, `'instructor'`, or
`'student'`. `classroom.ready` resolves when the client library and event handlers
are ready, not when a server connection or classroom join has succeeded.
`start()` and `join()` return promises that resolve after queuing their requests;
the role changes when the server acknowledges the request.

For an existing UI manager, use {@link StateManager#enableLearning} to create and
own the classroom controller:

```javascript
const classroom = controls.enableLearning({
    mode: 'hidden',
    session: 'biochem-bonding-class-a'
});
```

Calling `enableLearning` again replaces the previous connection;
`controls.enableLearning(false)` disables classroom support. Destroying the UI
manager or viewer also releases its classroom connection, controls, and listeners.
Applications that manage this lifecycle themselves can construct
`new $3Dmol.LearningSession(controls, options)` and must call
`classroom.destroy()` when removing it. Use one controller per viewer.

`classroom.leave()` leaves the classroom, ending it when called by its instructor.
In hidden mode, explicitly leaving also stops automatic joining. Call
`controls.enableLearning(options)` again to resume waiting for an instructor.
When the instructor ends a session, other hidden-mode viewers return to waiting
for the next activation of their configured session.

The Socket.IO browser client is loaded when classroom support is enabled. If
your page already supplies a compatible Socket.IO client, it is reused. Set
`LearningConfig.clientURL` to a self-hosted client script when a content security
policy or offline asset setup requires it. That script must provide a Socket.IO
4.x client; the session server remains a separate deployment.

### Python notebooks and exported HTML

The [py3Dmol constructor](https://github.com/3dmol/3Dmol.js/blob/master/py3Dmol/README.md#editing-trajectory-and-classroom-controls)
also enables these controls. Use `learning=True` for the classroom button, or
`learning={"mode": "hidden", "session": "my-unique-classroom"}` for hidden
activation. The same dictionary accepts `server`, `hotkey`, and `clientURL`.
Independent `edit=True` and `trajectory=True` options enable the other viewer
controls. py3Dmol loads the optional UI bundle automatically, and these options
work in both notebook displays and exported HTML. See its README for executable
examples and bundle configuration.

## Shared scenes

The instructor's edit-menu changes are shared with connected students: model
loading and removal, selections and representations, residue labels, molecular
surfaces, volumetric maps, distance measurements, background, projection, and trajectory
settings. Students can rotate and zoom locally, but joining closes their edit
panel and disables editing and trajectory controls until they leave the session.
An instructor camera change overrides their local camera view.

Use **Volumetric data** to load an OpenDX, Gaussian cube, CCP4/MRC, or VASP map
from a local file or URL. Compressed `.gz` files are supported. Choose an
isosurface and adjust its contour level, color, opacity, or visibility; the parsed
grid is reused when changing the contour. Direct volume rendering uses WebGL 2.
URL downloads require the remote server to permit cross-origin requests.

**Surfaces** provides molecular surface calculation and material controls.
Annotation shapes supplied in shared URLs are also restored and synchronized.
Students receive the resulting geometry and map data, so they do not need the
instructor's local files or repeat surface calculations.

Distance measurements are opt-in under **Measurements**; see the
[viewer controls](tutorial-url.html#distance-measurements) for creating and
removing them. Each measurement has color and label-visibility controls. The
instructor's measurements, appearance changes, and removals are shared with
students, whose editing controls remain disabled while connected.

Changes are grouped briefly before transmission, and only changed scene sections
are sent. Camera updates are separate from molecular data. Joining or reconnecting
retrieves the current scene; ordinary appearance changes preserve unchanged
models and geometry on student viewers. Local files used in a shared scene are
sent to the session server and connected viewers.

Closing the instructor's tab or window ends the session. Once the server detects
the disconnect, students' previous control settings are restored and the session name becomes
available again. Losing the instructor's connection also ends the session;
students can reconnect while the instructor remains connected.

Sessions also expire after **48 hours without activity**. Valid participant actions,
including instructor updates, queries, joins, and reconnections, reset this idle
period. The server checks expiry when a session is accessed and runs a cleanup
pass every minute. Expiry disconnects viewers from the session and makes its name
available again. Restarting the server also clears all sessions because they are
held in memory.

<iframe width="560" height="315" src="https://www.youtube.com/embed/90UhGlzLcdc" frameborder="0" allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>

## Self-Hosting

The easiest way to use the 3Dmol.js active learning environment is through the hosted viewer at [https://3dmol.org/viewer.html](https://3dmol.org/viewer.html).  However, it may be desirable to host your own instance for your class.  A light-weight Flask webserver can quickly be brought up on a Ubuntu Linux machine with a few commands:


```{@lang bash}
#install dependencies
apt install npm python3-pip python3-venv git
#get latest version from git
git clone https://github.com/3dmol/3Dmol.js.git
#build with npm
cd 3Dmol.js
python3 -m venv .venv
. .venv/bin/activate
python3 -m pip install -r learning-environment/requirements.txt
npm install
#run a standalone flask server
cd learning-environment
#can optionally specify a port with -p <PORT>, default is 5000
python3 server.py --host 0.0.0.0
```

The commands above are a setup outline. Install the repository's supported Node.js 24 and npm 10 or newer before building; a distribution's default `npm` package may not meet those requirements. The server requires Python 3.10 or newer and uses native threads with `simple-websocket` and Gunicorn, following the [Flask-SocketIO threaded deployment](https://flask-socketio.readthedocs.io/en/latest/deployment.html#gunicorn-web-server). Eventlet is no longer required. When upgrading an existing server, update its virtual environment from `learning-environment/requirements.txt` and restart it.

The command above listens on all interfaces at port 5000, e.g. `http://HOSTNAME:5000/`. Without `--host`, it binds only to `127.0.0.1`, suitable for a local reverse proxy. The viewer is at `http://HOSTNAME:5000/static/viewer.html`.
Your network must allow external connections to the specified port.
The server runs one Gunicorn worker with 100 threads; use `--threads` to adjust capacity. Keep one worker because session state is held in process memory. See the [session-server README](https://github.com/3dmol/3Dmol.js/blob/master/learning-environment/README.md) for the equivalent Gunicorn command.
Do <b>not</b> run this command as superuser to run the server on the privileged default webserver port 80 without consulting with a network security specialist.

The standalone viewer connects to the session server at the page's origin, preserving its hostname, protocol, and port. This applies equally to `localhost` and `127.0.0.1`; neither address selects port 5000 automatically. Embedded classrooms instead default to `https://3dmol.csb.pitt.edu`; point their `data-learning-server` or `LearningConfig.server` at your server. To serve the standalone viewer and session server separately, set `window.viewerSessionServer` before the viewer initializes, for example in a script in the viewer page:

```javascript
window.viewerSessionServer = 'http://localhost:6000';
```

Use the exact address and port on which the session server listens. In particular, `localhost` can resolve to IPv6, while the default server bind is IPv4 `127.0.0.1`. The browser also treats these hostnames as different origins. Remove the override, or set it to `window.location.origin`, when proxying through the webserver that serves the viewer.

Deploy the updated session server with the rebuilt UI bundle when enabling hidden
embeds. Automatic student discovery requires the server's session-watch events;
an older server can accept ordinary joins but will not announce instructor
activation to waiting viewers. Allow the embed origin to connect through your
server/proxy CORS policy, and use HTTPS for the session server when embedding on
an HTTPS page.

If port 5000 is occupied, start the server on another port, such as `python3 server.py -p 5068`, and open `http://127.0.0.1:5068/static/viewer.html`. This direct setup needs no Apache configuration. On macOS, a `403` response with a `Server: AirTunes` header indicates that the request reached Apple's service rather than the session server; adding CORS headers in Apache will not fix that response.

The `learning-environment/static` symlink points to the repository root. Files placed in that root can be referenced as relative paths within the `url=` component.  For example: `http://localhost:5000/static/viewer.html?url=molecule.sdf&style=stick`


If Apache serves the viewer, proxy `/socket.io/` to the session server and keep the default page-origin connection. Both HTTP polling and WebSocket upgrades must reach the same backend. Apache 2.4.47 and later can handle both through `mod_proxy` and `mod_proxy_http`; for a backend started with `-p 5068`, add this to the virtual host and reload Apache (see [Apache's protocol-upgrade documentation](https://httpd.apache.org/docs/2.4/mod/mod_proxy.html#protoupgrade)):

```apache.conf
ProxyPass        /socket.io/ http://127.0.0.1:5068/socket.io/ upgrade=websocket
ProxyPassReverse /socket.io/ http://127.0.0.1:5068/socket.io/
```

Use one proxy configuration matching your backend port. Older Apache setups can instead use the following rewrite-based example for a server on port 5000 (see the [Socket.IO proxy discussion](https://stackoverflow.com/questions/36472920/apache-proxy-configuration-for-socket-io-project-not-in-root)):


```apache.conf
RewriteEngine On
RewriteCond %{REQUEST_URI}  ^/socket.io            [NC]
RewriteCond %{QUERY_STRING} transport=websocket    [NC]
RewriteRule /(.*)           ws://localhost:5000/$1 [P,L]

ProxyPass        /socket.io http://localhost:5000/socket.io
ProxyPassReverse /socket.io http://localhost:5000/socket.io
```

The [flask-SocketIO documentation](https://flask-socketio.readthedocs.io/en/latest/) provides more information for running in other setups such as gunicorn, uwsgi, and nginx.

