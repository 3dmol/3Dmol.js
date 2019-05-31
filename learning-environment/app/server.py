from flask import Flask, render_template
from flask_socketio import SocketIO, send

app = Flask(__name__)
app.config['SECRET_KEY'] = 'secret!'
socketio = SocketIO(app)

@socketio.on('message')
def handleMessage(msg):
	print('Message: ' + msg)
    send(msg, broadcast=True)    
# @app.route('/')
# def function():
#     return app.send_static_file('viewer.html')


if __name__ == '__main__':
    socketio.run(app)