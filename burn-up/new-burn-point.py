import datetime

def proc():
    # f = open("../rel-todo-gtk3")
    f = open("../menu-items-to-fix")
    lines = f.readlines()
    add_to_pre_sum = True
    n_above = 0
    n_below = 0
    for line in lines:
        if "--- Done ---" in line:
            add_to_pre_sum = False
        arr = line.split()
        try:
            idx = arr.index(":::")
            if idx < len(arr) - 1:
                estimated_time_for_task_string = arr[idx+1]
                etft = float(estimated_time_for_task_string)
                if add_to_pre_sum:
                    n_above += etft
                else:
                    n_below += etft
        except ValueError as e:
            pass

    return (n_above, n_below, n_above + n_below)


(n_above, n_below, n_total) = proc()


start = datetime.datetime(2021, 9, 29)
now   = datetime.datetime.now()

delta_s = int(now.strftime("%s")) - int(start.strftime("%s"))
delta_d = delta_s/(60*60*24)
print(delta_d, n_above, n_below, n_total)
