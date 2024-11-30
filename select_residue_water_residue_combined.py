from pymol import cmd
import math

def select_residue_water_residue_combined(domainA, domainB):
    '''
    Identifies and lists residues from two domains that interact via bridging water molecules,
    considering both polar contacts and hydrogen bonds between residues and water molecules.
    The script differentiates between polar contacts and hydrogen bonds based on geometric criteria
    and lists distances and total distances in the results.

    Enhancements:
    - Excludes water molecules from residue selections.
    - Calculates and lists distances between residues and water molecules.
    - Calculates and lists the total distance from residue A to residue B via the water molecule.
    - Differentiates between polar contacts and hydrogen bonds in the visualization by using different colors.
    - Displays these interactions, including distance labels.

    Parameters:
    - domainA: The selection name or selection string for the first domain.
    - domainB: The selection name or selection string for the second domain.
    '''

    # Parameters for interactions
    max_total_distance = 7.0  # Maximum total distance for interactions
    max_residue_water_distance = max_total_distance / 2.0  # Approximate max distance between residue and water

    # Parameters for hydrogen bonds
    hbond_distance_cutoff = 3.5  # Maximum distance between donor and acceptor atoms for hydrogen bonds
    hbond_angle_cutoff = 120     # Minimum angle (degrees) for hydrogen bonds

    # Select all water molecules in the structure
    cmd.select('all_waters', 'resn HOH')

    # Identify water molecules within max_residue_water_distance of both domains
    cmd.select('water_near_A', 'all_waters within %.1f of (%s)' % (max_residue_water_distance, domainA))
    cmd.select('water_near_B', 'all_waters within %.1f of (%s)' % (max_residue_water_distance, domainB))
    cmd.select('water_bridge', 'water_near_A and water_near_B')

    # Debugging outputs
    num_waters_A = cmd.count_atoms('water_near_A and name O')
    num_waters_B = cmd.count_atoms('water_near_B and name O')
    num_bridging_waters = cmd.count_atoms('water_bridge and name O')
    print('Number of water molecules near domainA:', num_waters_A)
    print('Number of water molecules near domainB:', num_waters_B)
    print('Number of bridging water molecules:', num_bridging_waters)

    # Get list of bridging water molecules
    water_indices = cmd.index('water_bridge and name O and alt ""')

    if not water_indices:
        print('No bridging water molecules found between the specified domains within the specified distance.')
        cmd.delete('all_waters')
        cmd.delete('water_near_A')
        cmd.delete('water_near_B')
        cmd.delete('water_bridge')
        return

    interaction_list = []

    # Create selections for visualization
    cmd.select('interacting_residues_domain1', 'none')
    cmd.select('interacting_residues_domain2', 'none')
    cmd.select('interacting_waters', 'none')

    # Interaction IDs for unique naming
    hbond_id = 0
    polar_contact_id = 0

    for (water_obj, water_idx) in water_indices:
        water_sele = '(%s and index %d and name O and alt "")' % (water_obj, water_idx)

        # Get water coordinates
        water_coords = cmd.get_atom_coords(water_sele)

        # Check for hydrogen bonds between water and domainA residues
        hbonds_A = find_hydrogen_bonds(domainA, water_sele, hbond_distance_cutoff, hbond_angle_cutoff)

        # Check for hydrogen bonds between water and domainB residues
        hbonds_B = find_hydrogen_bonds(domainB, water_sele, hbond_distance_cutoff, hbond_angle_cutoff)

        # Check for polar contacts (excluding hydrogen bonds) between water and domain residues
        polar_contacts_A = []
        polar_contacts_B = []

        if not hbonds_A:
            polar_contacts_A = find_polar_contacts(domainA, water_sele, max_residue_water_distance)
        if not hbonds_B:
            polar_contacts_B = find_polar_contacts(domainB, water_sele, max_residue_water_distance)

        # Determine interaction type based on available data
        if hbonds_A and hbonds_B:
            # Water forms hydrogen bonds with both domainA and domainB residues
            # Process hydrogen bond interactions
            hbond_id = process_interactions(hbonds_A, hbonds_B, water_obj, water_idx,
                                            hbond_id, 'hbond', interaction_list,
                                            domainA, domainB, water_sele, water_coords)
        elif (hbonds_A and polar_contacts_B) or (polar_contacts_A and hbonds_B):
            # Mixed interactions (hydrogen bond with one domain and polar contact with the other)
            # Process as polar contacts
            contacts_A = hbonds_A or polar_contacts_A
            contacts_B = hbonds_B or polar_contacts_B
            polar_contact_id = process_interactions(contacts_A, contacts_B, water_obj, water_idx,
                                                    polar_contact_id, 'polar_contact', interaction_list,
                                                    domainA, domainB, water_sele, water_coords)
        elif polar_contacts_A and polar_contacts_B:
            # Water forms polar contacts with both domainA and domainB residues
            # Process polar contact interactions
            polar_contact_id = process_interactions(polar_contacts_A, polar_contacts_B, water_obj, water_idx,
                                                    polar_contact_id, 'polar_contact', interaction_list,
                                                    domainA, domainB, water_sele, water_coords)

    # Visualize the interactions
    if interaction_list:
        # Show and color residues from domain1
        cmd.show('sticks', 'interacting_residues_domain1')
        cmd.color('yellow', 'interacting_residues_domain1')

        # Show and color residues from domain2
        cmd.show('sticks', 'interacting_residues_domain2')
        cmd.color('magenta', 'interacting_residues_domain2')

        # Show and color water molecules
        cmd.show('spheres', 'interacting_waters')
        cmd.color('cyan', 'interacting_waters')
        cmd.set('sphere_scale', 0.25, 'interacting_waters')  # Set sphere size to 0.25

        # Set visualization properties for hydrogen bonds
        cmd.hide('labels', 'hbond_*')
        cmd.set('dash_color', 'green', 'hbond_*')
        cmd.set('dash_width', 2.0)
        cmd.show('dashes', 'hbond_*')
        cmd.label('hbond_*', '"%.2f Å" % (distance)')

        # Set visualization properties for polar contacts
        cmd.hide('labels', 'polar_contact_*')
        cmd.set('dash_color', 'orange', 'polar_contact_*')
        cmd.set('dash_width', 1.5)
        cmd.show('dashes', 'polar_contact_*')
        cmd.label('polar_contact_*', '"%.2f Å" % (distance)')
    else:
        print('No residue interactions found via bridging water molecules within the specified criteria.')

    # Output the total interactions and their details
    total_interactions = len(interaction_list)
    print('Total interactions: %d' % total_interactions)

    # List interactions with distances
    for interaction in interaction_list:
        print('Residue %s interacts with residue %s via water id %d (%s). Distances: %.2f Å (resA-water), %.2f Å (water-resB), Total distance: %.2f Å' %
              (interaction['resA'], interaction['resB'], interaction['water'], interaction['type'],
               interaction['distA'], interaction['distB'], interaction['total_distance']))

    # Clean up temporary selections
    cmd.delete('all_waters')
    cmd.delete('water_near_A')
    cmd.delete('water_near_B')
    cmd.delete('water_bridge')

def find_hydrogen_bonds(domain_sele, water_sele, distance_cutoff, angle_cutoff):
    '''
    Finds hydrogen bonds between a domain and a water molecule.
    Returns a list of atom pairs involved in hydrogen bonds.
    '''
    hbonds1 = cmd.find_pairs('(%s and (name N,O) and not resn HOH)' % domain_sele, water_sele,
                             mode=1, cutoff=distance_cutoff, angle=angle_cutoff)
    hbonds2 = cmd.find_pairs(water_sele, '(%s and (name N,O) and not resn HOH)' % domain_sele,
                             mode=1, cutoff=distance_cutoff, angle=angle_cutoff)
    return hbonds1 + hbonds2

def find_polar_contacts(domain_sele, water_sele, distance_cutoff):
    '''
    Finds polar contacts between a domain and a water molecule within a distance cutoff.
    Returns a list of atom pairs involved in the polar contacts.
    '''
    polar_atoms = cmd.get_model('(%s and (name N*,O*) and alt "" and not resn HOH)' % domain_sele).atom
    water_atom = cmd.get_model('%s and name O and alt ""' % water_sele).atom[0]

    polar_contacts = []
    for atom in polar_atoms:
        dist = distance_between_coords(atom.coord, water_atom.coord)
        if dist <= distance_cutoff:
            polar_contacts.append(((atom.model, atom.index), (water_atom.model, water_atom.index), dist))
    return polar_contacts

def process_interactions(contacts_A, contacts_B, water_obj, water_idx, interaction_id, interaction_type,
                         interaction_list, domainA, domainB, water_sele, water_coords):
    '''
    Processes interactions and updates selections and interaction lists.
    Calculates distances and total distances for each interaction.
    '''
    # Get residues involved and calculate distances
    resA_set = set()
    resA_distances = {}
    for contact in contacts_A:
        atom_idx = contact[0][1]
        atom = cmd.get_model('index %d' % atom_idx).atom[0]
        if atom.resn == 'HOH':
            continue  # Skip water molecules
        resA_info = (atom.model, atom.chain, atom.resi, atom.resn)
        resA_set.add(resA_info)
        # Calculate distance to water
        distA = distance_between_coords(atom.coord, water_coords)
        resA_distances[resA_info] = (atom_idx, distA)

    resB_set = set()
    resB_distances = {}
    for contact in contacts_B:
        atom_idx = contact[0][1]
        atom = cmd.get_model('index %d' % atom_idx).atom[0]
        if atom.resn == 'HOH':
            continue  # Skip water molecules
        resB_info = (atom.model, atom.chain, atom.resi, atom.resn)
        resB_set.add(resB_info)
        # Calculate distance to water
        distB = distance_between_coords(atom.coord, water_coords)
        resB_distances[resB_info] = (atom_idx, distB)

    if not resA_set or not resB_set:
        return interaction_id  # No valid interactions

    interaction_id += 1  # Increment interaction counter

    for resA_info in resA_set:
        for resB_info in resB_set:
            distA = resA_distances[resA_info][1]
            distB = resB_distances[resB_info][1]
            total_distance = distA + distB
            if total_distance <= 7.0:
                # Record the interaction
                interaction_list.append({
                    'resA': '%s/%s`%s/%s' % resA_info,
                    'resB': '%s/%s`%s/%s' % resB_info,
                    'water': water_idx,
                    'type': 'Hydrogen Bond' if interaction_type == 'hbond' else 'Polar Contact',
                    'distA': distA,
                    'distB': distB,
                    'total_distance': total_distance
                })

                # Add residues to the interacting_residues_domain1 and domain2 selections
                resA_residue_sele = 'model %s and chain %s and resi %s' % (resA_info[0], resA_info[1], resA_info[2])
                resB_residue_sele = 'model %s and chain %s and resi %s' % (resB_info[0], resB_info[1], resB_info[2])

                cmd.select('interacting_residues_domain1', 'interacting_residues_domain1 or (%s)' % resA_residue_sele)
                cmd.select('interacting_residues_domain2', 'interacting_residues_domain2 or (%s)' % resB_residue_sele)
                # Add water molecule to the interacting_waters selection
                cmd.select('interacting_waters', 'interacting_waters or %s' % water_sele)

                # Create distance objects
                atomA_sele = 'index %d' % resA_distances[resA_info][0]
                atomB_sele = 'index %d' % resB_distances[resB_info][0]
                if interaction_type == 'hbond':
                    cmd.distance('hbond_%d_A' % interaction_id, atomA_sele, water_sele)
                    cmd.distance('hbond_%d_B' % interaction_id, water_sele, atomB_sele)
                else:
                    cmd.distance('polar_contact_%d_A' % interaction_id, atomA_sele, water_sele)
                    cmd.distance('polar_contact_%d_B' % interaction_id, water_sele, atomB_sele)

    return interaction_id

def distance_between_coords(coord1, coord2):
    '''
    Calculates the Euclidean distance between two coordinates.
    '''
    return math.sqrt((coord1[0] - coord2[0]) ** 2 +
                     (coord1[1] - coord2[1]) ** 2 +
                     (coord1[2] - coord2[2]) ** 2)

# Extend PyMOL commands
cmd.extend('select_residue_water_residue_combined', select_residue_water_residue_combined)