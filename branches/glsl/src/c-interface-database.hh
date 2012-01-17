/* src/c-interface-database.hh
 * 
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

/*  ----------------------------------------------------------------------- */
/*             Database                                                     */
/*  ----------------------------------------------------------------------- */
/*! \brief set the host name of the database */
set_db_host(const char *host);
/*! \brief set the user name of the database */
set_db_username(const char *username);
/*! \brief the password name of the database */
set_db_password(const char *password);

int connect_database();

#endif // USE_MYSQL_DATABASE
