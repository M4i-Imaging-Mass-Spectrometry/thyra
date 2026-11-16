"""Pictures of a converted store, from the command line.

``thyra ion-image`` writes one ion image as a PNG and ``thyra report``
writes a one-file HTML overview. Both read a store that already exists and
change nothing in it.

Nothing is imported here: the commands import what they need when they
run, so ``thyra --version`` and the conversion command pay for none of it.

The modules, in the order a command uses them:

- :mod:`.store` opens the store and reads single pieces of it;
- :mod:`.ion_image` turns an m/z value into a window of columns and the
  per-pixel sums into a raster;
- :mod:`.render` turns a raster into pixels with numpy and Pillow;
- :mod:`.report` lays the pieces out as HTML;
- :mod:`.cli` holds the two commands.
"""
