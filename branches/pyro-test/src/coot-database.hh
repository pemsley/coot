
// #define USE_MYSQL_DATABASE

#ifdef USE_MYSQL_DATABASE
#include <mysql/mysql.h>
#include <mysql/errmsg.h>

void db_finish_up();
int db_query_insert(const std::string &insert_string);

#endif

