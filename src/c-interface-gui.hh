

#include <gtk/gtk.h>
#include <string>

// this might need a better name

// get string for column 0 (which are strings)
std::string get_active_label_in_combobox(GtkComboBox *combobox);

void pepflips_by_difference_map_dialog();

void set_transient_for_main_window(GtkWidget *dialog);

class ProgressBarPopUp {
	GtkWindow* window;
	GtkProgressBar* progress_bar;

	public:
	ProgressBarPopUp(const char* title, const char* description);
    ProgressBarPopUp(const ProgressBarPopUp&) = delete;
    ProgressBarPopUp(ProgressBarPopUp&&);
    
	void pulse();
	void set_fraction(float frac);
	~ProgressBarPopUp();
};