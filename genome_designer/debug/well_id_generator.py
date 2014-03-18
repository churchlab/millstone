"""
Generator object for standard 96-well plate.
"""


class WellIdGenerator(object):
    """Generates 96-plate well ids from A1 ... A12, ..., H1 ... H12

    Also returns plate number if requested.
    """

    LETTER_TRANSITION_TABLE = {
            'A': 'B',
            'B': 'C',
            'C': 'D',
            'D': 'E',
            'E': 'F',
            'F': 'G',
            'G': 'H',
            'H': 'A',
    }


    def __init__(self, include_plate=False):
        self.letter = 'A'
        self.number = 1
        self.include_plate = include_plate
        if include_plate:
            self.plate = 1


    def __iter__(self):
        return self


    def next(self):
        # Create the current return value.
        current_id = self.letter + "%02d" % (self.number,)

        # Bump the state.
        if self.number == 12:
            self.letter = self.LETTER_TRANSITION_TABLE[self.letter]
            # If we are back to A, bump the plate number.
            if self.include_plate and self.letter == 'A':
                self.plate += 1

        if self.number == 12:
            self.number = 1
        else:
            self.number += 1

        # Return current.
        if self.include_plate:
            return (self.plate, current_id)
        return current_id
