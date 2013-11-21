# GBVTD

## Ground Based Velocity Track Display

This tarball contains three directories:

  * cappi/ : 'cappi' creator using intebilinear interpolation
  * findcenter/ : 'findcenter' simplex search code
  * vtd/ : 'VD' GBVTD wind retrieval code

To build, run 

     $ make -f Makefile.(compiler)

in each directory. 
You will need an input file for each
radar volume you which to process for VD, while findcenter supports multiple
volumes in the same input file. The input files are commented.

Notes about the Linux version:
The Linux port reads and writes files in cedric ASCII format, with the
header readable by 'Bigendian' machines (ie Solaris, HP). This is done so that
existing applications (ie cedric, grid2ps) can still read the files. It cannot
currently read or write binary cedric files,but this is planned for the future.

Good luck on the build and resulting analysis!

## Contributing to GBVTD

* Check out the latest master to make sure the feature hasn't been implemented or the bug hasn't been fixed yet
* Check out the [issue tracker](http://github.com/mmbell/gbvtd/issues) to make sure someone already hasn't requested it and/or contributed it
* Fork the project
* Start a feature/bugfix branch
* Commit and push until you are happy with your contribution, then open a Pull Request.
* Make sure to add tests for the feature/bugfix. This is important so I don't break it in a future version unintentionally.
* Please try not to mess with the Rakefile, version, or history. If you want to have your own version, or is otherwise necessary, that is fine, but please isolate it to its own commit so I can cherry-pick around it.

## Copyright

Copyright (c) 2013 Michael Bell, Wen-Chau Lee

See LICENSE for details.

