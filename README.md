## Updating Boost sources

I chose the set of submodules by checking out boost.math and:

1. Try to build. If it succeeds, stop.
2. Otherwise, a header file is missing. Search a full boost enlistment to determine which library includes that header (we'll use the name "foo" here), and add it to the set of submodules:
   1. `git submodule add foo https://github.com/boostorg/foo`
   2. `git submodule init foo`

Checkout all of the submodules at the tag for the appropriate Boost version:
```
git submodule for_each "git checkout boost-1.66.0"
```

## TODO
* Build, build, build!
* special_math.cpp needs the "new sources" procedure.
* plumb in the test
