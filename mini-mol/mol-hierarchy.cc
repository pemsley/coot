
	 for(unsigned int ifrag=0; ifrag<target_ca_coords.fragments.size(); ifrag++) {
	    for(int ires=target_ca_coords[ifrag].min_res_no(); ires<=target_ca_coords[ifrag].max_residue_number(); ires++) {
	       for (unsigned int iat=0; iat<target_ca_coords[ifrag][ires].atoms.size(); iat++) {
		  std::cout << " " << target_ca_coords[ifrag].fragment_id << " " << ires << " " << target_ca_coords[ifrag][ires]
			    << " " << target_ca_coords[ifrag][ires][iat].name
			    << " " << target_ca_coords[ifrag][ires][iat].pos.format() << std::endl;
	       }
	    }
	 }
