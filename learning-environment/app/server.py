from flask import Flask, request
from flask_socketio import SocketIO, send, emit, join_room, leave_room, close_room 
import json

# for socketio
import eventlet
eventlet.monkey_patch()

app = Flask(__name__)
app.config['SECRET_KEY'] = 'secret!'
socketio =  SocketIO(app, async_mode='eventlet')

sessions = {} # indexed by session name


@socketio.on('message')
def handleMessage(msg):
    print('Message: ' + msg)
    send(msg, broadcast=True)

@socketio.on('check session name event')
def handleCheckName(json):
    name = json['name']
    print("checking ",name,name in sessions)
    if name in sessions:
        emit('check session name response', "exists")
    else:
        emit('check session name response', "free")    

@socketio.on('create session event')
def handleCreateSession(json):
    name = json['name']
    if name in sessions:
        print("Invalid name",name)
        emit('create session response', False)
    else: 
        emit('create session response', True)
        print("Session Created",name)
        sessions[name] = {'cnt': 0,
                          'sid': request.sid, # id of the controller
                          'state': json['state'],
                          'view': json['view'],
                          'query_updates': {}, # each clients individual selections
                          'query_active' : False
                             }        
        join_room(name)
        emit('session count', 0, room=request.sid)

@socketio.on('join session event')
def handleJoinSession(json):
    name = json['name']
    print("Join",name)
    if name not in sessions:
        emit('join session response', "0")
    else: 
        emit('join session response', "1")
        join_room(name)
        try:
            sessions[name]['cnt'] += 1
        except:
            emit('error: restart connection')
        emit('session count', sessions[name]['cnt'], room=sessions[name]['sid'])
        #send session state to connecting client
        if 'view' in sessions[name]:
            emit('viewer view change response', sessions[name]['view'], room=request.sid)
        if 'state' in sessions[name]:
            emit('viewer state change response', sessions[name]['state'], room=request.sid)
        print("Joined Session: " + str(json))
        
        if sessions[name]['query_active']:
            emit('query start response', sessions[name]['query_active'], room=request.sid)  

@socketio.on('delete session event')
def handleDeleteSession(json):
    name = json['name']
    if name not in sessions:
        return
    if request.sid != sessions[name]['sid']:
        return #only the controller can send updates    
    emit('leave session response', room=name)
    close_room(name)
    print("Session Deleted:" + str(json))
    del sessions[name]
    emit('delete session response')



@socketio.on('leave session event')
def handleLeaveSession(json):
    name = json['name']
    leave_room(name)
    if name not in sessions:
        return
    sessions[name]['cnt'] -= 1
    
    updates = sessions[name]['query_updates']
    if request.sid in updates:
        del updates[request.sid]
        emit('query update response', len(updates), room=name)   

    emit('session count', sessions[name]['cnt'], room=sessions[name]['sid'])
    print("Session Left: " + str(json))
    emit('leave session response')


@socketio.on('viewer view change event')
def handleViewerChange(json):
    name = json['name']
    viewstate = json['view']
    if name not in sessions:
        return
    if request.sid != sessions[name]['sid']:
        print("Incorrect sid for viewer change")
        return #only the controller can send updates
    sessions[name]['view'] = viewstate
    emit('viewer view change response', viewstate, room=name)
    
@socketio.on('viewer state change event')
def handleStateChange(json):
    name = json['name']
    state = json['state']
    print("State change",name)
    if name not in sessions:
        return
    if request.sid != sessions[name]['sid']:
        print("Incorrect sid for state change")
        return #only the controller can send updates    
    sessions[name]['state'] = state
    emit('viewer state change response', state, room=name)    

@socketio.on('query start')
def handleQueryStart(json):
    name = json['name']
    print("query start",name)
    if name not in sessions:
        return
    if request.sid != sessions[name]['sid']:
        return #only the controller can send updates    
    sessions[name]['query_updates'] = {}
    sessions[name]['query_active'] = json['kind']

    emit('query start response', json['kind'], room=name)    
    
@socketio.on('query end')
def handleQueryEnd(json):
    name = json['name']
    if name not in sessions:
        return
    if request.sid != sessions[name]['sid']:
        return #only the controller can send updates    
    sessions[name]['query_updates'] = {}
    sessions[name]['query_active'] = False
    emit('query end response', '', room=name) 
    
@socketio.on('query update')
def handleQueryUpdate(json):
    '''A single update from a client, respond with number unique responses'''
    name = json['name']
    sel = json['selected']
    if name not in sessions:
        return
    # process update     
    updates = sessions[name]['query_updates']    
    updates[request.sid] = sel
    
    emit('query update response', len(updates), room=sessions[name]['sid'])   
    
@socketio.on('query fetch')
def handleQueryFetch(json):
    '''Return the merged results'''
    name = json['name']
    if name not in sessions:
        return
    # process update     
    updates = sessions[name]['query_updates']
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
    socketio.run(app)
