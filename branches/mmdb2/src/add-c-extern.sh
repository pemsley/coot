

if [ -z "$1" ] ; then
  echo usage $0 filename
  exit 1
fi

# ok to proceed with adding

mv "$1" "$1".tmp

echo "#ifndef BEGIN_C_DECLS"                 > "$1".tmp.a 
echo "#ifdef __cplusplus"                   >> "$1".tmp.a 
echo "#define BEGIN_C_DECLS extern \"C\" {" >> "$1".tmp.a 
echo "#define END_C_DECLS }"                >> "$1".tmp.a 
echo "#else"                                >> "$1".tmp.a 
echo "#define BEGIN_C_DECLS"                >> "$1".tmp.a 
echo "#define END_C_DECLS"                  >> "$1".tmp.a 
echo "#endif"                               >> "$1".tmp.a 
echo "#endif"                               >> "$1".tmp.a
echo ""              >> "$1".tmp.a
echo "BEGIN_C_DECLS" >> "$1".tmp.a 
echo ""              >> "$1".tmp.a

echo ""                                                             > "$1".tmp.b
# echo "#if (((GTK_MAJOR_VERSION == 2) && (GTK_MINOR_VERSION > 5)) || GTK_MAJOR_VERSION > 2)" >> "$1".tmp.b
echo "#if (GTK_MAJOR_VERSION > 1)"                                 >> "$1".tmp.b
echo "char* coot_revision(void);"                                  >> "$1".tmp.b
echo "GtkWidget* create_aboutdialog (void);"                       >> "$1".tmp.b
echo "GtkWidget* create_coords_filechooserdialog1(void);"          >> "$1".tmp.b
echo "GtkWidget* create_dataset_filechooserdialog1(void);"         >> "$1".tmp.b
echo "GtkWidget* create_map_name_filechooserdialog1(void);"        >> "$1".tmp.b
echo "GtkWidget* create_phs_coordinates_filechooserdialog1(void);" >> "$1".tmp.b
echo "GtkWidget* create_save_coords_filechooserdialog1(void);"     >> "$1".tmp.b
echo "GtkWidget* create_cif_dictionary_filechooserdialog1(void);"  >> "$1".tmp.b
echo "GtkWidget* create_run_script_filechooserdialog1(void);"      >> "$1".tmp.b
echo "GtkWidget* create_save_symmetry_coords_filechooserdialog1(void);" >> "$1".tmp.b
echo "GtkWidget* create_save_state_filechooserdialog1(void);"      >> "$1".tmp.b
echo "GtkWidget* create_screendump_filechooserdialog1(void);"      >> "$1".tmp.b
echo "GdkPixbuf* create_pixbuf(const gchar *filename);"            >> "$1".tmp.b
echo "gchar* find_pixmap_file (const gchar     *filename);"        >> "$1".tmp.b
echo "GtkWidget* create_model_toolbar_menu (void);"                >> "$1".tmp.b
echo "GtkWidget         *create_residue_editor_select_monomer_type_dialog();" >> "$1".tmp.b
echo "GtkWidget *wrapped_create_residue_editor_select_monomer_type_dialog();" >> "$1".tmp.b
echo "GtkWidget* create_restraints_editor_dialog (void);"                     >> "$1".tmp.b
echo "GtkWidget* create_residue_editor_select_monomer_type_dialog (void);"    >> "$1".tmp.b
echo "GtkWidget* create_save_restraint_chooserdialog (void);"                 >> "$1".tmp.b
# this probably shouldnt be here, but I dont want to add this function to
# gtk1 version
echo "GtkWidget* create_run_refmac_nolabels_help_dialog(void);"    >> "$1".tmp.b
echo "GtkWidget* create_run_refmac_file_help_dialog(void);"        >> "$1".tmp.b
echo "GtkWidget* create_run_refmac_sad_help_dialog(void);"         >> "$1".tmp.b
echo "GtkWidget* create_run_refmac_mtz_filechooserdialog(void);"   >> "$1".tmp.b
echo "GtkWidget* create_coot_references_dialog(void);"             >> "$1".tmp.b
echo "GtkWidget* create_fast_ss_search_dialog(void);"              >> "$1".tmp.b
echo "GtkWidget* create_baton_build_params_dialog(void);"          >> "$1".tmp.b
echo "GtkWidget* create_pisa_interfaces_dialog(void);"             >> "$1".tmp.b
echo "GtkWidget* create_scheme_window(void);" 	                   >> "$1".tmp.b
# echo "void on_accession_code_ok_button_clicked(GtkButton *button, gpointer user_data);" >> "$1".tmp.b
echo "#endif /* GTK_MAJOR_VERSION */ "                             >> "$1".tmp.b
echo ""            >> "$1".tmp.b
echo "END_C_DECLS" >> "$1".tmp.b
echo ""            >> "$1".tmp.b

cat "$1".tmp.a  "$1".tmp "$1".tmp.b > "$1"

rm  "$1".tmp.a  "$1".tmp "$1".tmp.b

