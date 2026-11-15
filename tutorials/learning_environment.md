<script src="../build/3Dmol-min.js"></script> 

The [hosted 3Dmol.js viewer](tutorial-url.html) can serve as an active learning environment of engaging and querying students when teaching properties of 3D molecular structures.  The instructor creates a session using the share button in the upper right corner.  When students connect to this session, they see the same molecular view as the instructor.  Although the students can independently manipulate the view, the instructor can override their view.

The instructor can pose questions that the students then answer by clicking on atoms.  The accumulated and anonymous answers can then be displayed in the instructor's view.

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
the disconnect, students regain editing controls and the session name becomes
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

The viewer connects to the session server at the page's origin, preserving its hostname, protocol, and port. This applies equally to `localhost` and `127.0.0.1`; neither address selects port 5000 automatically. To serve the viewer and session server separately, set `window.viewerSessionServer` before `initSessions()` runs, for example in a script in the viewer page:

```javascript
window.viewerSessionServer = 'http://localhost:6000';
```

Use the exact address and port on which the session server listens. In particular, `localhost` can resolve to IPv6, while the default server bind is IPv4 `127.0.0.1`. The browser also treats these hostnames as different origins. Remove the override, or set it to `window.location.origin`, when proxying through the webserver that serves the viewer.

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

