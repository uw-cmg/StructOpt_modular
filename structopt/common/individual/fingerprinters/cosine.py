from scipy.spatial.distance import cosine


def cosine(fingerprint1, fingerprint2):
    """Calculates the cosine distance / 2.0 between two fingerprints."""
    return cosine(fingerprint1, fingerprint1) / 2.0

