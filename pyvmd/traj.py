"""
Trajectory analysis utilities
"""
import logging

from .objects import Molecule

__all__ = ['Loader', 'LoadStatus']


LOGGER = logging.getLogger(__name__)


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
        self.molecule.frame = self._chunk_frame


class Loader(object):
    """
    Iteratively loads the trajectory files.
    """
    def __init__(self, molecule, traj_files, step=1, chunk=10):
        """
        @param molecule: Molecule used for loading the trajectory.
        @param traj_files: List of trajectory files
        @param step: Load every 'step'th frame from trajectory.
        @param chunk: Number of frames to load at once
        """
        assert isinstance(molecule, Molecule)
        assert isinstance(step, (int, long)) and step > 0, "step must be possitive integer"
        assert isinstance(chunk, (int, long)) and chunk > 0, "chunk must be possitive integer"
        self.molecule = molecule
        self.traj_files = traj_files
        self.step = step
        self.chunk = chunk
        self._callbacks = []

    def add_callback(self, callback):
        """
        Add callback to be called on every frame.
        """
        self._callbacks.append(callback)

    def add_collector(self, collector):
        """
        Adds collector to be run on every frame.
        """
        self._callbacks.append(collector.collect)

    def run(self):
        """
        Run the trajectory.
        """
        # Clear the molecule frames
        del self.molecule.frames[:]

        status = LoadStatus(self.molecule)
        for filename in self.traj_files:
            start = 0
            while True:
                # Load 'chunk' frames
                stop = start + self.step * self.chunk - 1
                LOGGER.debug('Loading %s from %d to %d, every %d', filename, start, stop, self.step)
                self.molecule.load(filename, start=start, stop=stop, step=self.step)
                loaded = len(self.molecule.frames)
                if not loaded:
                    # No frames were loaded
                    break

                # Call the callback
                status.next_chunk()
                for dummy in xrange(0, loaded):
                    status.next_frame()
                    LOGGER.debug('Analyzing frame %d', status.frame)
                    for callback in self._callbacks:
                        callback(status)

                # Prepare for next iteration - delete all frames
                del self.molecule.frames[:]
                if loaded < self.chunk:
                    # Nothing else to be loaded for this filename
                    break
                start += self.step * self.chunk
        LOGGER.info('Analyzed %s frames.', status.frame + 1)
