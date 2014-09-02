/* src/c-interface-database.cc
 * 
 * Copyright 2001, 2002, 2003, 2004, 2005, 2006, 2007 by The University of York
 * Copyright 2007 by The University of Oxford
 * 
 * Author: Paul Emsley
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA
 */

#ifdef USE_MYSQL_DATABASE

#include "c-interface-database.hh"

/*  ----------------------------------------------------------------------- */
/*             Database                                                     */
/*  ----------------------------------------------------------------------- */

/*! \brief set the host name of the database */
set_db_host(const char *host) {
   graphics_info_t::mysql_host = host;
}

/*! \brief set the user name of the database */
set_db_username(const char *username) {
   graphics_info_t::mysql_user = username;
}

/*! \brief the password name of the database */
set_db_password(const char *password) {
   graphics_info_t::mysql_passwd = password;
}

int connect_database() {

   std::string host = "localhost";
   std::string user = "cootuser";
   std::string passwd = "password";

   const char *db = "cootsessions";
   unsigned int port = 0;
   const char *unix_socket = 0; 
   unsigned long client_flag = 0;

   // graphics_info_t::mysql = new MYSQL; // not needed.
   graphics_info_t::mysql = mysql_init(graphics_info_t::mysql);

   if (graphics_info_t::mysql) { 
      graphics_info_t::mysql = mysql_real_connect(graphics_info_t::mysql,
						  host.c_str(), user.c_str(), passwd.c_str(),
						  db, port, unix_socket, client_flag);

      // We need to do a test here to see if that connection worked.
      // 
      if (graphics_info_t::mysql == 0) {

	 std::cout << "INFO:: Can't connect to database." << std::endl;

      } else { 

	 // We only want to do this if this person doesn't exist
	 // already...
	 std::string query("insert into user (userid, username) values");
	 std::pair<std::string, std::string> p = coot::get_userid_name_pair();
	 graphics_info_t::db_userid_username = p;
	 query += " ('";
	 query += p.first;
	 query += "','";
	 query += p.second;
	 query += "');";
	 unsigned long length = query.length();
      
	 int v = mysql_real_query(graphics_info_t::mysql,
				  query.c_str(), length);

	 if (v != 0) {
	    if (v == CR_COMMANDS_OUT_OF_SYNC)
	       std::cout << "WARNING:: MYSQL Commands executed in an"
			 << " improper order" << std::endl;
	    if (v == CR_SERVER_GONE_ERROR) 
	       std::cout << "WARNING:: MYSQL Server gone!"
			 << std::endl;
	    if (v == CR_SERVER_LOST) 
	       std::cout << "WARNING:: MYSQL Server lost during query!"
			 << std::endl;
	    if (v == CR_UNKNOWN_ERROR) 
	       std::cout << "WARNING:: MYSQL Server transaction had "
			 << "an uknown error!" << std::endl;
	 }
	 // std::cout << "start: mysql_real_query returned " << v << std::endl;

	 // set the session-id
	 time_t *timep = new time_t;
	 *timep = time(0); 
	 char *ct = ctime(timep);
	 delete timep;
	 std::string sct(ct);
	 // now remove the irritating carriage return at the end of the
	 // ctime string
	 std::string::size_type icarr = sct.find("\n");
	 if (icarr != std::string::npos) {
	    sct = sct.substr(0,icarr);
	 } 
	 std::string sessionid("session-");
	 sessionid += sct;
	 graphics_info_t::sessionid += sessionid;
      }

   } else {
      std::cout << "WARNING:: can't init mysql database structure"
		<< std::endl;
   }
}


#endif // USE_MYSQL_DATABASE
