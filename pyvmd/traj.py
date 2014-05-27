"""
Trajectory analysis utilities
"""
import logging

from Molecule import Molecule


class LoadStatus(object):
    """
    Holds information about loading.

    @ivar molecule: Molecule object
    @ivar frame: Total count of current frame
    """
    def __init__(self, molecule):
        self.molecule = molecule
        # Total frame count
        self.frame = -1
        # Frame number of the currently loaded frame
        self._chunk_frame = -1

    def __repr__(self):
        return '<%s: %d>' % (type(self).__name__, self.frame)

    def __str__(self):
        return 'Frame %d' % self.frame

    def next_chunk(self):
        """
        New chunk is loaded.
        """
        self._chunk_frame = -1

    def next_frame(self):
        """
        Move to the next frame.
        """
        self.frame += 1
        self._chunk_frame += 1
        self.molecule.setFrame(self._chunk_frame)


class Loader(object):
    """
    Iteratively loads the trajectory files.
    """
    def __init__(self, parm_file, traj_files, step=10):
        """
        @param parm_file: Structure file
        @param traj_files: List of trajectory files
        @param step: Number of frames to load at once
        """
        assert step > 0, "Number of steps must be possitive"
        self.molecule = Molecule()
        self.molecule.load(parm_file)
        self.traj_files = traj_files
        self.step = step
        self._callbacks = []

    def add_callback(self, callback):
        """
        Add callback to be called on every frame.
        """
        self._callbacks.append(callback)

    def run(self):
        """
        Run the trajectory.
        """
        status = LoadStatus(self.molecule)
        for filename in self.traj_files:
            start = 0
            while True:
                # Load 'step' frames
                stop = start + self.step - 1
                logging.debug('Loading %s from %d to %d', filename, start, stop)
                self.molecule.load(filename, first=start, last=stop)
                loaded = self.molecule.numFrames()
                if not loaded:
                    # No frames were loaded
                    break

                # Call the callback
                status.next_chunk()
                for dummy in xrange(0, loaded):
                    status.next_frame()
                    for callback in self._callbacks:
                        callback(status)

                # Prepare for next iteration
                self.molecule.delFrame()
                if loaded < self.step:
                    # Nothing else to be loaded for this filename
                    break
                start += self.step
        logging.info('Analyzed %s frames.', status.frame + 1)
