#!/usr/bin/python
#registers the use of 3Dmol.js in a database
#this script gets pinged everytime 3dmol.js is loaded,
#this allows us to report usage to our funding agencies
import os,MySQLdb

host = os.environ["REMOTE_ADDR"]
domain = host
if "HTTP_REFERER" in os.environ:
   domain = os.environ["HTTP_REFERER"]

conn = MySQLdb.connect (host = "localhost",user = "tracker",db="3dmoltrack")
cursor = conn.cursor()
cursor.execute("INSERT INTO accessed (host,domain) VALUES(%s,%s)",(host,domain[:111]))
conn.commit()
cursor.close()
conn.close()

print "Content-Type: text/html"     # HTML is following
#print "Access-Control-Allow-Origin: *" # allow cross-site scripting - disabled since enabled on server
print                               # blank line, end of headers
