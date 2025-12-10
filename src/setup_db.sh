#check if port is open
port=1738
ss -ltnu | grep -qE "[:.]$port\b" && echo "Port $port is IN USE" || echo "Port $port is FREE"

#setup db in a new folder
mkdir ./sock
initdb -D ./deca_db
pg_ctl -D ./deca_db -l ./deca_db/logfile start -o "-p $port -h 127.0.0.1 -k ~/WorkForaging/Academia/Nicole/deca/sock"
pg_isready -h 127.0.0.1 -p $port
createdb -h 127.0.0.1 -p $port deca_db 
createuser -h 127.0.0.1 -p $port -s postgres
psql -h 127.0.0.1 -p $port -d postgres -c "ALTER DATABASE deca_db OWNER TO postgres"
