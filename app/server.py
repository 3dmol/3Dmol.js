from flask import Flask, render_template

app = Flask(__name__)

# @app.route("/")
# def home():
#     return render_template("index.html")

@app.route("/viewer.html")
def viewer():
    return render_template("viewer.html")

if __name__ == "__main__":
    app.run(debug=True)