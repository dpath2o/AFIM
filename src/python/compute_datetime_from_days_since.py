from datetime import datetime, timedelta 
day_amount = 729024.5                                     # Sample day data construction
start_datetime = datetime(1,1,1,0,0,0)                 # Sample start datetime construction
print(start_datetime)
print(start_datetime + timedelta(days = day_amount))      # Generating and printing the end datetime
# 2007-10-19 00:00:00