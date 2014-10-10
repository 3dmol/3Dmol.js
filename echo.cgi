#!/usr/local/bin/python
#fetchs a file - workaround for cross site scripting prevention
import cgi,requests,sys


form = cgi.FieldStorage()
url = form["url"].value

response = requests.get(url)
print "Content-Type: text/text"
print ""

sys.stdout.write(response.text)

