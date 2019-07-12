Dependencies
####
pip3 install flask flask_socketio eventlet

TODOs
#####
JSHint the session code
Implement save restore of state for labels, shapes, and surfaces.
Improve efficiency of state updates
 - make sure updates are only sent if something has changed
 - only send changed part
 - could record and playback a history (this is most efficient)
Add ability to select residues
  -Query Atoms button becomes a pulldown button that can be either Query Atoms or Residues
Change style of clicked label on hover
Refactor code to be general and module
  - interface shoudl be separate from sessions information
  - create a Session class
  - use Promise API to handle event/response
Implement garbage collecting of stale sessions (24h)
Implement status reporting/concise logging
Implement password based authentication of controller (NOT in viewer.html)
Implement ability to create session not as controller (initializes state, but does not control)
Uses Session API to implement environment within py3Dmol
 -this will use the password based authentication to allow instructor to take control
