# Interface reference

```@autodocs
Modules = [ MolecularMinimumDistances ]
Order = [ :function, :struct ]
Filter = (f) -> (nameof(f) in (
    :minimum_distances,
    :minimum_distances!,
    :MinimumDistance,
    :AllPairs,
    :SelfPairs,
    :CrossPairs,
))
```