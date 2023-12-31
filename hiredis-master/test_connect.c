#include <stdlib.h>
#include "hiredis.h"

int main (int argc, char **argv){
    const char *hostname = (argc > 1) ? argv[1] : "127.0.0.1";
    int port = (argc > 2)? atoi(argv[2]) :6379; //6379 is the default port for redis...
    unsigned isunix =  0;
    redisContext *c = redisConnect (hostname, port);
    struct timeval timeout = {1, 50000}; //1.5 seconds
    
    if (isunix){
        c = redisConnectUnixWithTimeout (hostname, timeout);
        //providing time to connect with the server
    }
    else{
        c = redisConnectWithTimeout (hostname, port, timeout);
    }

    if (c == NULL || c -> err) {
        if (c){
            printf("Connection error : %s\n", c->errstr);
            //freeing the memory is connection not established successfully
            redisFree(c);
        }
        else {
            printf("Connection error : cannot allocate redis context\n");
        }

        exit(1);
    }
printf ("Connected to redis\n");
}