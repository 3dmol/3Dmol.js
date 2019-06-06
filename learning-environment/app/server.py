from flask import Flask
from flask_socketio import SocketIO, send, emit, join_room, leave_room, close_room
import json

# for socketio
import eventlet
eventlet.monkey_patch()

app = Flask(__name__)
app.config['SECRET_KEY'] = 'secret!'
socketio =  SocketIO(app, async_mode='eventlet')

session_dictionary = {} #stores session_name : session_password
session_count = {}

def check_session_present(session_json):
    if session_json['session_name'] in session_dictionary: #searches in O(1)
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
        session = session_json['session_name']
        password = session_json['session_password']
        session_dictionary[session] = password     
        join_room(session)
        session_count[session] = 0
        emit('session count', session_count[session], room=session)


@socketio.on('delete session event')
def handleDeleteSession(session_json):
    session = session_json['session_name']
    emit('leave session response', room=session)
    close_room(session)
    print("Session Deleted:" + str(session_json))
    session_dictionary.pop(session)
    session_count.pop(session)
    emit('delete session response')

@socketio.on('join session event')
def handleJoinSession(session_json):
    if (check_session_present(session_json) == False):
        emit('join session response', "0")
    else: 
        session = session_json['session_name']
        emit('join session response', "1")
        join_room(session)
        try:
            session_count[session] += 1
        except:
            emit('error: restart connection')
        emit('session count', session_count[session], room=session)
        print("Joined Session: " + str(session_json))

@socketio.on('leave session event')
def handleLeaveSession(session_json):
    session = session_json['session_name']
    leave_room(session)
    session_count[session] -= 1
    emit('session count', session_count[session], room=session)
    print("Session Left: " + str(session_json))
    emit('leave session response')

if __name__ == '__main__':
    socketio.run(app)