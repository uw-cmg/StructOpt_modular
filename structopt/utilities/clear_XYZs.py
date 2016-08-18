import os
import shutil

def clear_XYZs(n, generation, path):
    """Depending on the value in the post_processing dictionary, clear old
    XYZ files to save space. Run this every generation to clean up XYZ files

    Parameters
    ----------
    n : int
        if n < 0 : Only generation up to current generation - n are kept
        if n > 0 : Every n generation is kept

    generation : int
        The generation currently on. This function always keeps this generation file.
    path : string
        The /path/to/folder that holds the XYZs folder
    """

    assert(type(n) is int)

    # If nothing is getting removed (on the first generation)
    if generation == 0:
        return

    # Keeping the last n generations
    if n < 0 and self.generation > -n:
        path = os.path.join(path, 'XYZs/generation{}'.format(generation + n))
    # Keeping every n generation
    elif n > 0 and generation > 1 and generation % n != 1:
        path = os.path.join(path, 'XYZs/generation{}'.format(generation - 1))

    if path is not None:
        shutil.rmtree(path)

    return
