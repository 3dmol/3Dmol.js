#!/usr/bin/env python3
#fetchs a file - workaround for cross site scripting prevention
import cgi,requests,sys


form = cgi.FieldStorage()
url = form["url"].value

response = requests.get(url)
sys.stdout.write("Content-Type: text/text\n\n")
sys.stdout.write(response.text)

