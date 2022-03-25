# Advanced usage

### Advanced, non-allocating call

If this function will be called from within a loop where the coordinates of the molecules involved change, a faster and non-allocating call can be used, by exposing the interface of `CellListMap`. For example:

```julia-repl
julia> function iterate_lists(nsteps, list, water, protein, box, cl, aux)
           inds(i) = mol_index(i,3)
           for i in 1:nsteps
               # water and protein coordinates, and box, could change here
               cl = UpdateCellList!(water, protein, box, cl, aux)
               minimum_distances!(inds, list, box, cl)
               # list was updated
           end
        end


```


