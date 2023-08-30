<script src="../build/3Dmol-min.js"></script> 

The [hosted 3Dmol.js viewer](tutorial-url.html) can serve as an active learning environment of engaging and querying students when teaching properties of 3D molecular structures.  The instructor creates a session using the share button in the upper right corner.  When students connect to this session, they see the same molecular view as the instructor.  Although the students can independently manipulate the view, the instructor can override their view.

The instructor can pose questions that the students then answer by clicking on atoms.  The accumulated and anonymous answers can then be displayed in the instructor's view.

<iframe width="560" height="315" src="https://www.youtube.com/embed/90UhGlzLcdc" frameborder="0" allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>

## Self-Hosting

The easiest way to use the 3Dmol.js active learning environment is through the hosted viewer at [https://3dmol.csb.pitt.edu].  However, it may be desirable to host your own instance for your class.  A light-weight Flask webserver can quickly be brought up on a Ubuntu Linux machine with a few commands:


```{@lang bash}
#install dependencies
apt install npm python3-pip git
pip3 install flask flask_socketio eventlet
#get latest version from git
git clone https://github.com/3dmol/3Dmol.js.git
#build with npm
cd 3Dmol.js
npm install
#run a standalone flask server
cd learning-environment
#can optionally specify a port with -p <PORT>, default is 5000
python3 server.py
```

By default, this will launch a server on port 5000, e.g. `http://HOSTNAME:5000/`.  The viewer is at `http://HOSTNAME:5000/static/viewer.html`.
Your network must allow external connections to the specified port.
Do <b>not</b> run this command as superuser to run the server on the privileged default webserver port 80 without consulting with a network security specialist.

Note that files put in the `3dmol` directory can be referenced as relative paths within the `url=` component.  For example: `http://localhost:5000/static/viewer.html?url=molecule.sdf&style=stick`


If you already have a webserver running, you can redirect the websocket requests used by 3Dmol.js to the server running on port 5000.  For example, for apache the following is added to the configuration (see [here](https://stackoverflow.com/questions/36472920/apache-proxy-configuration-for-socket-io-project-not-in-root)</a>):


<pre>
<font color='black'>RewriteEngine On
RewriteCond %{REQUEST_URI}  ^/socket.io            [NC]
RewriteCond %{QUERY_STRING} transport=websocket    [NC]
RewriteRule /(.*)           ws://localhost:5000/$1 [P,L]

ProxyPass        /socket.io http://localhost:5000/socket.io
ProxyPassReverse /socket.io http://localhost:5000/socket.io
<font>
</pre>

The [flask-SocketIO documentation](https://flask-socketio.readthedocs.io/en/latest/) provides more information for running in other setups such as gunicorn, uwsgi, and nginx.

