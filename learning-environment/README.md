# Active learning session server

This Flask/Socket.IO service synchronizes classroom viewers and collects atom
selections. The browser client is `sessions.js`; `server.py` keeps sessions in
memory, so restarting the server discards them. Ordinary embedded viewers do not
need this service.

For setup, viewer URLs, and proxy configuration, use the maintained
[active learning tutorial](../tutorials/learning_environment.md). Build the main
repository first, then run `python3 server.py` from this directory; `-p` selects
the port (default 5000). The `static` symlink exposes the repository root, and `/`
redirects to `/static/viewer.html`.

## Dependencies

Install these in a Python virtual environment:

pip3 install flask flask_socketio eventlet

## Session-server backlog

These are historical development ideas, not a list of supported features or a
release plan. Verify each against the current server, client, and UI tests before
starting work.

Lint the session code as part of the current JavaScript development workflow
Provide viewergrid functionality with the ability to query molecules
Add support for volumetric data in viewer
Add support for uploading your own files (nontrivial)
Implement save restore of state for labels, shapes, and surfaces.
Improve efficiency of state updates
 - make sure updates are only sent if something has changed
 - only send changed part
 - could record and playback a history (this is most efficient)
Add ability to select residues
  -Query Atoms button becomes a pulldown button that can be either Query Atoms or Residues
Refactor code to be general and module
  - interface should be separate from sessions information
  - create a Session class
  - use Promise API to handle event/response
Implement garbage collecting of stale sessions (24h)
Implement status reporting/concise logging
Implement password based authentication of controller (NOT in viewer.html)
Implement ability to create session not as controller (initializes state, but does not control)
Uses Session API to implement environment within py3Dmol
 -this will use the password based authentication to allow instructor to take control
