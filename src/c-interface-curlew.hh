
void curlew_install_extension_file(const std::string &file_name, const std::string &checksum,
                                   GtkWidget *install_button, GtkWidget *uninstall_button);

// return uninstall status (true for done)

bool curlew_uninstall_extension_file(const std::string &file_name);

