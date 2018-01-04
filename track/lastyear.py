#!/usr/bin/python
#extract statistics for last year
import os,MySQLdb


conn = MySQLdb.connect (host = "localhost",user = "tracker",db="3dmoltrack")
cursor = conn.cursor()
cursor.execute("SELECT DISTINCT(domain) FROM `accessed`  WHERE `time` >= DATE_SUB(NOW(),INTERVAL 1 YEAR)")

for row in cursor:
    print row

