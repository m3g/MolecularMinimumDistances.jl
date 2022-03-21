#
# Functions used to reduce the lists in parallel calculations
#

"""
```
keep_best_list!(list1,list2)
```

Update `list1` considering the list of minimum distances given in `list2`.

Internal function or structure - interface may change. 

"""
function keep_best_list!(list1, list2)
    for i in eachindex(list1)
        if list2[i].d < list1[i].d
            list1[i] = list2[i]
        end
    end
end

"""

```
reduce_list!(list, list_threaded)
```

Update the final `list` of minimum-distances given the threaded list `list_threaded`.

Internal function or structure - interface may change. 

"""
function reduce_list!(list, list_threaded)
    list .= list_threaded[1]
    for it in 2:length(list_threaded)
        keep_best_list!(list, list_threaded[it])
    end
    return list
end

"""

```
reduce_list_pair!(lists, lists_threaded)
```

Update the final tuple `lists`, of minimum-distances of the two molecule sets given the threaded lists `lists_threaded`.

Internal function or structure - interface may change. 

"""
function reduce_list_pair!(lists, lists_threaded)
    lists[1] .= lists_threaded[1][1]
    lists[2] .= lists_threaded[1][2]
    for it in 2:length(lists_threaded)
        keep_best_list!(lists[1], lists_threaded[it][1])
        keep_best_list!(lists[2], lists_threaded[it][2])
    end
    return lists
end
