import datetime

def proc():
    # f = open("../rel-todo-gtk3")
    fn = "../menu-items-to-fix"
    fn = "../rel-todo-gtk4"
    f = open(fn)
    lines = f.readlines()
    # add_to_other_items = False
    add_to_other_items = True
    add_to_done_other_items = False
    n_other_items = 0
    n_done_other_items = 0
    for line in lines:
        if "--- Items ---" in line:
            add_to_other_items = True
            add_to_done_other_items = False
        if "--- Done Items ---" in line:
            add_to_other_items = False
            add_to_done_other_items = True
        arr = line.split()
        try:
            idx = arr.index(":::")
            if idx < len(arr) - 1:
                estimated_time_for_task_string = arr[idx+1]
                etft = float(estimated_time_for_task_string)
                if add_to_other_items:
                    n_other_items += etft
                if add_to_done_other_items:
                    n_done_other_items += etft
        except ValueError as e:
            pass

    return (n_other_items, n_done_other_items)


(n_other_items, n_done_other_items) = proc()


start = datetime.datetime(2022, 8, 6)
now   = datetime.datetime.now()

delta_s = int(now.strftime("%s")) - int(start.strftime("%s"))
delta_d = delta_s/(60*60*24)
# print(delta_d, n_above, n_below, n_total)
# date/time n-items-not-yet done n-items-done
print("{:.3f} {} {}".format(delta_d, n_other_items, n_done_other_items))
