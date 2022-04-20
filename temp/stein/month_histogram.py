import mysql.connector
import numpy as np
import matplotlib.pyplot as plt
import math
from matplotlib.patches import Wedge


# NB Set up SSH tunnel first (forward remote port 3306 to local port 3309):  
# ssh -p 899 -L3309:127.0.0.1:3306 mmsteam@81.169.221.160 -N &

def main():
    
    conn = mysql.connector.connect(host="127.0.0.1", port=3309, user="mmsteam", db="MMS")
    cursor = conn.cursor()
    
    SQL = "SELECT month(DateStart),count(*) from MMS2up where DateUpdated > '2017-04-24' group by month(DateStart);"
    cursor.execute(SQL)
    result = cursor.fetchall()

    month = []
    freq = []
    for column in result:
        month.append(column[0])
        freq.append(column[1])


# --- plotting -- 

    t12=[1,2,3,4,5,6,7,8,9,10,11,12]
    my_xticks = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
    plt.xticks(t12, my_xticks)


    plt.ylabel('Number of events')
    plt.bar(month,freq)
    plt.show()

main()
