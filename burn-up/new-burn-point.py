import datetime

def proc():
    # f = open("../rel-todo-gtk3")
    f = open("../menu-items-to-fix")
    lines = f.readlines()
    add_to_menu_items = False
    add_to_other_items = False
    add_to_done_menu_items = False
    add_to_done_other_items = False
    n_menu_items  = 0
    n_other_items = 0
    n_done_menu_items  = 0
    n_done_other_items = 0
    for line in lines:
        if "--- Menu Items ---" in line:
            add_to_menu_items = True
            add_to_done_menu_items = False
            add_to_other_items = False
            add_to_done_other_items = False
        if "--- Done Menu Items ---" in line:
            add_to_menu_items = False
            add_to_done_menu_items = True
            add_to_other_items = False
            add_to_done_other_items = False
        if "--- Other Items ---" in line:
            add_to_menu_items = False
            add_to_done_menu_items = False
            add_to_other_items = True
            add_to_done_other_items = False
        if "--- Done Other Items ---" in line:
            add_to_menu_items = False
            add_to_done_menu_items = False
            add_to_other_items = False
            add_to_done_other_items = True
        arr = line.split()
        try:
            idx = arr.index(":::")
            if idx < len(arr) - 1:
                estimated_time_for_task_string = arr[idx+1]
                etft = float(estimated_time_for_task_string)
                if add_to_menu_items:
                    n_menu_items += etft
                if add_to_done_menu_items:
                    n_done_menu_items += etft
                if add_to_other_items:
                    n_other_items += etft
                if add_to_done_other_items:
                    n_done_other_items += etft
        except ValueError as e:
            pass

    return (n_menu_items, n_done_menu_items, n_other_items, n_done_other_items)


(n_menu_items, n_done_menu_items, n_other_items, n_done_other_items) = proc()


start = datetime.datetime(2021, 9, 29)
now   = datetime.datetime.now()

delta_s = int(now.strftime("%s")) - int(start.strftime("%s"))
delta_d = delta_s/(60*60*24)
# print(delta_d, n_above, n_below, n_total)
print("{:.3f} {} {:.2f} {} {}".format(delta_d, n_menu_items, n_done_menu_items, n_other_items, n_done_other_items))
