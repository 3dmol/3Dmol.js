from flask import Flask
from flask_socketio import SocketIO, send, emit, join_room, leave_room
import json

# for socketio
import eventlet
eventlet.monkey_patch()

app = Flask(__name__)
app.config['SECRET_KEY'] = 'secret!'
socketio =  SocketIO(app, async_mode='eventlet')

session_json_list = [] #stores session json
session_name_list = [] #stores session name

def check_session_present(session_json):
    if session_json['session_name'] in session_name_list:
        return True
    else:        
        return False

@socketio.on('message')
def handleMessage(msg):
    print('Message: ' + msg)
    send(msg, broadcast=True)

@socketio.on('check session name event')
def handleCheckName(session_json):
    if (check_session_present(session_json) == True):
        emit('check session name response', "0")
    else:
        emit('check session name response', "1")    

@socketio.on('create session event')
def handleCreateSession(session_json):
    if (check_session_present(session_json) == True):
        print("Invalid name")
        emit('create session response', "0")
    else: 
        emit('create session response', "1")
        print("Session Created: " + str(session_json))
        session_json_list.append(session_json)
        session = session_json['session_name']
        session_name_list.append(session)
        join_room(session)
   
    
@socketio.on('join session event')
def handleJoinSession(session_json):
    session = session_json['session_name']
    join_room(session)
    print("received my event: " + str(session_json))

if __name__ == '__main__':
    socketio.run(app)