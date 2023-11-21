import numpy as np

from __functions import *


def check_sunflowers_rapeseed_RotCnstr(vals_):
    crop_ = crop_names_dict['sunflowers']
    n_years = len(vals_)
    # before and after sunflowers there should not be rapeseed
    for year_ in range(n_years):
        if year_ >= 1 and year_ < n_years - 2:
            if crop_ == crop_names_dict['sunflowers']:
                values_suflowers = np.where(np.array(vals_) == crop_, 1, 0)
                values_rapeseed = np.where(np.array(vals_) == crop_names_dict['rapeseed'], 1, 0)
                values = values_suflowers + values_rapeseed
                # this is the same function just for the original input data
                if values[year_] > 0.5 and (values[year_ - 1] > 0.5 or values[year_ + 1] > 0.5):
                    # constraint violated -> return 1
                    return 1
        elif year_ == n_years - 2:
            values_suflowers = np.where(np.array(vals_) == crop_, 1, 0)
            values_rapeseed = np.where(np.array(vals_) == crop_names_dict['rapeseed'], 1, 0)
            values = values_suflowers + values_rapeseed
            # this is the same function just for the original input data

            if values[year_] > 0.5 and (values[year_ - 1] > 0.5):
                # constraint violated -> return 1
                return 1
    return 0


def check_rapeseed_RotCnstr(vals_):
    crop_ = crop_names_dict['rapeseed']
    values = np.where(np.array(vals_) == crop_, 1, 0)
    n_years = len(vals_)
    for year_ in range(n_years):
        if year_ < n_years-2:
            # Enforce rapeseed constraint only if rapeseed in year AND in any of the preceding
            # two years
            # this is the same function just for the original input data
            if values[year_] > 0.5 and (values[year_ + 1] > 0.5 or values[year_ + 2] > 0.5):
                # constraint violated -> return 1
                return 1
        elif year_ == n_years - 2:
            # this is the same function just for the original input data
            if values[year_] > 0.5 and values[year_ + 1] > 0.5:
                # constraint violated -> return 1
                return 1
    return 0


def check_potato_RotCnstr(vals_):
    crop_ =  crop_names_dict['potato']
    for year_ in range(len(vals_)):
        if year_ >= 3:
            if crop_ == crop_names_dict['potato']:
                # Enforce rapeseed constraint only if rapeseed in year AND in any of the preceding
                # two years
                values = np.where(np.array(vals_) == crop_, 1, 0)
                # this is the same function just for the original input data
                if values[year_] > 0.5 and any(values[year_ - n] > 0.5 for n in range(1, 4)):
                    # constraint violated -> return 1
                    return 1
    return 0


def check_beet_RotCnstr(vals_):
    crop_ = crop_names_dict['beets']
    for year_ in range(len(vals_)):
        if year_ >= 3:
            if crop_ == crop_names_dict['beets']:
                # Enforce rapeseed constraint only if rapeseed in year AND in any of the preceding
                # two years
                values = np.where(np.array(vals_) == crop_, 1, 0)
                # this is the same function just for the original input data
                if values[year_] > 0.5 and any(values[year_ - n] > 0.5 for n in range(1, 4)):
                    # constraint violated -> return 1
                    return 1
    return 0


def check_legume_RotCnstr(vals_):
    crop_ = crop_names_dict['legumes']
    for year_ in range(len(vals_)):
        if year_ >= 1:
            if crop_ == crop_names_dict['legumes']:
                # Enforce rapeseed constraint only if rapeseed in year AND in any of the preceding
                # two years
                values = np.where(np.array(vals_) == crop_, 1, 0)
                if sum(values) > 3:
                    print('impossible to change')
                # this is the same function just for the original input data
                if values[year_] > 0.5 and any(values[year_ - n] > 0.5 for n in range(1, 2)):
                    # constraint violated -> return 1
                    return 1
    return 0


def longest_cereal_seq(vals_):
    crop_ = crop_names_dict['winter_cereals']
    if crop_ == crop_names_dict['winter_cereals']:
        values = np.where(np.array(vals_) == crop_, 1, 0)
        longest_seq, max_start_index = longest_sequence(values)
        return longest_seq


def longest_maize_seq(vals_):
    crop_ = crop_names_dict['maize']
    if crop_ == crop_names_dict['maize']:
        values = np.where(np.array(vals_) == crop_, 1, 0)
        longest_seq, max_start_index = longest_sequence(values)
        return longest_seq


def check_CropRotRules(historic_croptypes_dict):
    # here we check for each field if our predefined crop rotation rules were violated;
    # i need a function that returns a dict in the following structure:
    # constr_is_enforcable_dict {field_1: [1, 1, 1], field_2: [1, 0, 1],...}
    # for each field we check if each of our e.g. three constraints is enforced by the farmer; if yes -> if no -> 0;
    # this dict will be used with the lazyconstraints; if dict 0 -> don't enforce this constraint
    print('checking rules')
    print(historic_croptypes_dict)
    rapeseed_violation_list = []
    potato_violation_list = []
    beet_violation_list = []
    legume_violation_list = []
    sunflower_rp_violation_list = []
    cereals_longest_seq_list = []
    maize_longest_seq_list = []

    for key, value in historic_croptypes_dict.items():
        rapeseed_violation = check_rapeseed_RotCnstr(value)
        rapeseed_violation_list.append(rapeseed_violation)

        potato_violation = check_potato_RotCnstr(value)
        potato_violation_list.append(potato_violation)

        beet_violation = check_beet_RotCnstr(value)
        beet_violation_list.append(beet_violation)

        legume_violation = check_legume_RotCnstr(value)
        legume_violation_list.append(legume_violation)

        sunflower_rp_violation = check_sunflowers_rapeseed_RotCnstr(value)
        sunflower_rp_violation_list.append(sunflower_rp_violation)

        cereals_longest_seq_list.append(longest_cereal_seq(value))
        maize_longest_seq_list.append(longest_maize_seq(value))

    print(sum(rapeseed_violation_list), 'violations of the rapeseed constraint')
    print(sum(potato_violation_list), 'violations of the potato constraint')
    print(sum(beet_violation_list), 'violations of the beet constraint')
    print(sum(legume_violation_list), 'violations of the legumes constraint')
    print(sum(sunflower_rp_violation_list), 'violations of the sunflower-rapeseed constraint')

    violation_dict = {id_: (val1, val2, val3, val4, val5) for id_, val1, val2, val3, val4, val5 in zip(list(historic_croptypes_dict.keys()), rapeseed_violation_list, potato_violation_list,
                                                                beet_violation_list, legume_violation_list, sunflower_rp_violation_list)}
    longest_seq_dict = {id_: (val1, val2) for id_, val1, val2 in zip(list(historic_croptypes_dict.keys()),
                                                                     maize_longest_seq_list, cereals_longest_seq_list)}
    print('length of maize list', len(maize_longest_seq_list))
    return violation_dict, longest_seq_dict



