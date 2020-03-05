#!/usr/bin/env python3
#extract statistics for last year
import os,MySQLdb


conn = MySQLdb.connect (host = "localhost",user = "tracker",db="3dmoltrack")
cursor = conn.cursor()
cursor.execute("SELECT SUBSTRING_INDEX(SUBSTRING_INDEX(SUBSTRING_INDEX(SUBSTRING_INDEX(domain, '/', 3), '://', -1), '/', 1), '?', 1) FROM `accessed` WHERE `time` >= DATE_SUB(NOW(),INTERVAL 1 YEAR)")

for row in cursor:
    print(row[0])

