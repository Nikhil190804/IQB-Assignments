# Chou and Fasman parameters
alpha_propensities = {
    'E': 1.53, 'A': 1.45, 'L': 1.34, 'H': 1.24, 'M': 1.20, 'Q': 1.17, 'W': 1.14, 'V': 1.14,
    'F': 1.12, 'K': 1.07, 'I': 1.00, 'D': 0.98, 'T': 0.82, 'S': 0.79, 'R': 0.79, 'C': 0.77,
    'N': 0.73, 'Y': 0.61, 'P': 0.59, 'G': 0.53
}

beta_propensities = {
    'M': 1.67, 'V': 1.65, 'I': 1.60, 'C': 1.30, 'Y': 1.29, 'F': 1.28, 'Q': 1.23, 'L': 1.22,
    'T': 1.20, 'W': 1.19, 'A': 0.97, 'R': 0.90, 'G': 0.81, 'D': 0.80, 'K': 0.74, 'S': 0.72,
    'H': 0.71, 'N': 0.65, 'P': 0.62, 'E': 0.26
}

# Protein sequence
sequence = "MNASSEGESFAGSVQIPGGTTVLVELTPDIHICGICKQQFNNLDAFVAHKQSGCQLTGTSAAAPSTVQFVSEETVPATQTQTTTRTITSETQTITVSAPEFVFEHGYQTYLPTESNENQTATVISLPAKSRTKKPTTPPAQKRLNCCYPGCQFKTAYGMKDMERHLKIHTGDKPHKCEVCGKCFSRKDKLKTHMRCHTGVKPYKCKTCDYAAADSSSLNKHLRIHSDERPFKCQICPYASRNSSQLTVHLRSHTASELDDDVPKANCLSTESTDTPKAPVITLPSEAREQMATLGERTFNCCYPGCHFKTVHGMKDLDRHLRIHTGDKPHKCEFCDKCFSRKDNLTMHMRCHTSVKPHKCHLCDYAAVDSSSLKKHLRIHSDERPYKCQLCPYASRNSSQLTVHLRSHTGDTPFQCWLCSAKFKISSDLKRHMIVHSGEKPFKCEFCDVRCTMKANLKSHIRIKHTFKCLHCAFQGRDRADLLEHSRLHQADHPEKCPECSYSCSSAAALRVHSRVHCKDRPFKCDFCSFDTKRPSSLAKHVDKVHRDEAKTENRAPLGKEGLREGSSQHVAKIVTQRAFRCETCGASFVRDDSLRCHKKQHSDQSENKNSDLVTFPPESGASGQLSTLVSVGQLEAPLESQDL"

def calculate_propensity(window, propensity_dict):
    score = sum(propensity_dict.get(aa, 0) for aa in window)
    return score / len(window)

def chou_fasman(sequence):
    helix_regions = []
    beta_regions = []

    window_size_helix = 6  # Window size for helix propensity calculation
    window_size_strand = 5 # Window size for beta strand propensity calculation

    for i in range(len(sequence) - window_size_helix + 1):
        window_helix = sequence[i:i+window_size_helix]
        helix_score = calculate_propensity(window_helix, alpha_propensities)

        if helix_score > 1.03:  # Threshold for helix propensity
            helix_regions.append((i, i+window_size_helix))

    for i in range(len(sequence) - window_size_strand + 1):
        window_strand = sequence[i:i+window_size_strand]
        beta_score = calculate_propensity(window_strand, beta_propensities)

        if beta_score > 1.05:   # Threshold for beta strand propensity
            beta_regions.append((i, i+window_size_strand))

    return helix_regions, beta_regions

def resolve_conflicts(helix_regions, beta_regions, sequence):
    secondary_structure = [0] * len(sequence)

    for start, end in helix_regions:
        for i in range(start, end):
            if secondary_structure[i] == 0 or secondary_structure[i] == 1:
                secondary_structure[i] = 1
            elif secondary_structure[i] == 2:
                # Resolve conflict, choose the one with higher propensity score
                helix_window = sequence[start:end]
                beta_window = sequence[i:i+5]
                helix_score = calculate_propensity(helix_window, alpha_propensities)
                beta_score = calculate_propensity(beta_window, beta_propensities)
                if helix_score > beta_score:
                    for j in range(start, end):
                        secondary_structure[j] = 1
                else:
                    for j in range(i, i+5):
                        secondary_structure[j] = 2

    for start, end in beta_regions:
        for i in range(start, end):
            if secondary_structure[i] == 0 or secondary_structure[i] == 2:
                secondary_structure[i] = 2
            elif secondary_structure[i] == 1:
                # Resolve conflict, choose the one with higher propensity score
                helix_window = sequence[i:i+6]
                beta_window = sequence[start:end]
                helix_score = calculate_propensity(helix_window, alpha_propensities)
                beta_score = calculate_propensity(beta_window, beta_propensities)
                if beta_score > helix_score:
                    for j in range(start, end):
                        secondary_structure[j] = 2
                else:
                    for j in range(i, i+6):
                        secondary_structure[j] = 1

    return secondary_structure
def display_secondary_structure(sequence, secondary_structure):
    for i, ss in enumerate(secondary_structure):
        if i % 60 == 0:
            print()
        if ss == 1:
            print('H', end='')
        elif ss == 2:
            print('S', end='')
        else:
            print(' ', end='')
    print()
# Function to calculate average propensity
def avg_propensity(sequence, propensities):
    total = 0
    for residue in sequence:
        if residue in propensities:
            total += propensities[residue]
    return total / len(sequence)

# Function to predict secondary structure
def predict_structure(sequence):
    helix_windows = []
    strand_windows = []

    # Scan for valid windows
    window_size_helix = 6
    window_size_strand = 5
    for i in range(len(sequence) - window_size_helix + 1):
        window = sequence[i:i+window_size_helix]
        alpha_prop = avg_propensity(window, alpha_propensities)
        beta_prop = avg_propensity(window, beta_propensities)

        if alpha_prop > 1.03 and alpha_prop > beta_prop:
            helix_windows.append((i, i+window_size_helix))

    for i in range(len(sequence) - window_size_strand + 1):
        window = sequence[i:i+window_size_strand]
        alpha_prop = avg_propensity(window, alpha_propensities)
        beta_prop = avg_propensity(window, beta_propensities)

        if beta_prop > 1.05 and beta_prop > alpha_prop:
            strand_windows.append((i, i+window_size_strand))

    # Extend windows to the left and right
    extended_helix_windows = []
    extended_strand_windows = []

    for start, end in helix_windows:
        left = start
        right = end
        while left > 0 and avg_propensity(sequence[left-1:end], alpha_propensities) > 1.03:
            left -= 1
        while right < len(sequence) and avg_propensity(sequence[left:right+1], alpha_propensities) > 1.03:
            right += 1
        extended_helix_windows.append((left, right))

    for start, end in strand_windows:
        left = start
        right = end
        while left > 0 and avg_propensity(sequence[left-1:end], beta_propensities) > 1.05:
            left -= 1
        while right < len(sequence) and avg_propensity(sequence[left:right+1], beta_propensities) > 1.05:
            right += 1
        extended_strand_windows.append((left, right))

    # Create a list with 0, 1, or 2 for each position
    zero_one_list = [0] * len(sequence)

    for start, end in extended_helix_windows:
        for i in range(start, end):
            zero_one_list[i] = 1

    for start, end in extended_strand_windows:
        for i in range(start, end):
            if zero_one_list[i] == 1:
                zero_one_list[i] = 3  # Conflict
            else:
                zero_one_list[i] = 2

    # Resolve conflicts
    for i in range(len(zero_one_list)):
        if zero_one_list[i] == 3:
            window = sequence[max(0, i-3):i+4]
            alpha_prop = avg_propensity(window, alpha_propensities)
            beta_prop = avg_propensity(window, beta_propensities)
            if alpha_prop > beta_prop:
                zero_one_list[i] = 1
            else:
                zero_one_list[i] = 2

    # Print the secondary structure for helices
    print(" HELICAL IN NATURE\n")
    print("-" * 27)
    helix_regions = ""
    helix_pattern = ""
    for i, value in enumerate(zero_one_list):
        if value == 1:
            helix_regions += sequence[i]
            helix_pattern += "H"
        else:
            helix_regions += sequence[i] + ""
            helix_pattern += " "
        if (i + 1) % 80 == 0:
            print(helix_regions)
            print(helix_pattern)
            helix_regions = ""
            helix_pattern = ""
    if helix_regions:
        print(helix_regions)
        print(helix_pattern)
    print("\n" + "-" * 27 + "\n")

    # Print the secondary structure for strands
    print("BETA STRANDS IN NATURE\n")
    print("-" * 27)
    strand_regions = ""
    strand_pattern = ""
    for i, value in enumerate(zero_one_list):
        if value == 1:
            strand_regions += sequence[i]
            strand_pattern += "S"
        else:
            strand_regions += sequence[i] + ""
            strand_pattern += " "
        if (i + 1) % 80 == 0:
            print(strand_regions)
            print(strand_pattern)
            strand_regions = ""
            strand_pattern = ""
    if strand_regions:
        print(strand_regions)
        print(strand_pattern)
    print("\n" + "-" * 27 + "\n")

# Now let's call the function to predict the secondary structure
predict_structure(sequence)

helix_regions, beta_regions = chou_fasman(sequence)

secondary_structure = resolve_conflicts(helix_regions, beta_regions, sequence)
display_secondary_structure(sequence, secondary_structure)