
#include<string>

 extern "C" {

    GtkWidget *get_image_widget_for_comp_id(const std::string &comp_id, int imol, unsigned int image_size);
    GtkWidget *test_get_image_widget_for_comp_id(const std::string &comp_id);

 }
    
