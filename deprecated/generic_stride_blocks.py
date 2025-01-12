import numpy as np

def prop_identical_sliding_stride_np(s1, s2, w, stride):
    if w > len(s1):
        raise ValueError("Window size cannot be greater than sequence length")
    
    s1 = np.array(list(s1))
    s2 = np.array(list(s2))
    
    windows_s1 = np.array([s1[i:i+w] for i in range(0, len(s1) - w + 1, stride)])
    windows_s2 = np.array([s2[i:i+w] for i in range(0, len(s2) - w + 1, stride)])
    
    matches = np.sum(np.all(windows_s1 == windows_s2, axis=1))
    
    total_windows = len(windows_s1)
    return matches

# Example usage
s1 = "AATGTAGTGTCGTGCCGA"
s2 = "AATGTAGTGTAGTGCCGA"
w = 3
stride = 1
result = prop_identical_sliding_stride_np(s1, s2, w, stride)
print(f"Proportion of identical overlapping blocks: {result}")

