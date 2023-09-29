

#ifndef C_INTERFACE_GUI_HH
#define C_INTERFACE_GUI_HH
#include <gtk/gtk.h>
#include <string>
#include <memory>

// this might need a better name

// get string for column 0 (which are strings)
std::string get_active_label_in_combobox(GtkComboBox *combobox);

void pepflips_by_difference_map_dialog();

void set_transient_for_main_window(GtkWidget *dialog);

class ProgressBarPopUp {
	GtkWindow* window;
	GtkProgressBar* progress_bar;

	public:
	ProgressBarPopUp(const char* title, const char* description) noexcept;
    ProgressBarPopUp(const ProgressBarPopUp&) = delete;
    ProgressBarPopUp(ProgressBarPopUp&&) noexcept;
    
	void pulse() noexcept;
	void set_fraction(float frac) noexcept;
	void set_text(const char* text) noexcept;
	~ProgressBarPopUp();
};

/// A thread-safe wrapper around `ProgressBarPopUp`
class ProgressNotifier {
	std::shared_ptr<ProgressBarPopUp> progress_bar_popup;

	public:

	ProgressNotifier(std::shared_ptr<ProgressBarPopUp> popup) noexcept;
	/// This allows for reporting progress in a thread-safe way
	void update_progress(float frac) noexcept;
	/// This allows for reporting progress in a thread-safe way
	void pulse() noexcept;
	/// This allows for changing the label in a thread-safe way
	void set_text(const char* text) noexcept;

	~ProgressNotifier();

};


#endif // C_INTERFACE_GUI_HH