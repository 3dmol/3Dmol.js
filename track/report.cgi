#!/usr/bin/python
#registers the use of 3Dmol.js in a database
#this script gets pinged everytime 3dmol.js is loaded,
#this allows us to report usage to our funding agencies
import os,MySQLdb

host = os.environ["REMOTE_ADDR"];

conn = MySQLdb.connect (host = "localhost",user = "tracker",db="3dmoltrack")
cursor = conn.cursor()
cursor.execute("INSERT INTO accessed (host) VALUES(%s)",(host,))
conn.commit()
cursor.close()
conn.close()

print "Content-Type: text/html"     # HTML is following
print                               # blank line, end of headers
