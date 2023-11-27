import matplotlib.pyplot as plt
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


def check_crop_presuc(value, crop_t):
    # following crops of any crop type
    preceding_crop_list = []
    succeeding_crop_list = []
    v_seq = np.where(np.array(value) == crop_names_dict[crop_t], 1, 0)
    if sum(v_seq) > 0:
        for i, value_year in enumerate(v_seq):
            if value_year > 0:
                try:
                    if i-1 >= 0:
                        preceding = value[i - 1]
                        preceding_crop_list.append(crop_names_dict_reversed[preceding])
                except:
                    None

                try:
                    succeeding = value[i + 1]
                    succeeding_crop_list.append(crop_names_dict_reversed[succeeding])
                except:
                    None
    return preceding_crop_list, succeeding_crop_list


def get_rotations(historic_croptypes_dict):
    # this function extracts crop rotation characteristics such as the most frequent following crops
    final_df = pd.DataFrame(
        {'value': [], 'crop_t_from': [],
         'pre_suc': [], 'fid': []})
    for key_crop, value_crop in crop_names_dict.items():
        preceding_crop_sublist = []
        succeeding_crop_sublist = []
        ids_pre_list = []
        ids_suc_list = []
        for key, value in historic_croptypes_dict.items():
            if key == 0.0:
                continue
            rape_pre, rape_suc = check_crop_presuc(value, key_crop)
            preceding_crop_sublist.append(rape_pre)
            succeeding_crop_sublist.append(rape_suc)
            ids_pre_list.append(list(np.repeat(key, repeats=len(rape_pre))))
            ids_suc_list.append(list(np.repeat(key, repeats=len(rape_suc))))

        preceding_crop_sublist = sum(preceding_crop_sublist, [])
        preced_df = pd.DataFrame({'value': preceding_crop_sublist, 'crop_t_from': np.repeat(key_crop, repeats=len(preceding_crop_sublist)),
                                  'pre_suc': np.repeat("pre", repeats=len(preceding_crop_sublist)),
                                  'fid': sum(ids_pre_list, [])})

        final_df = final_df.append(preced_df)
        succeeding_crop_sublist = sum(succeeding_crop_sublist, [])
        succseed_df = pd.DataFrame({'value':succeeding_crop_sublist , 'crop_t_from': np.repeat(key_crop, repeats=len(succeeding_crop_sublist)),
                                  'pre_suc': np.repeat("suc", repeats=len(succeeding_crop_sublist)),
                                    'fid': sum(ids_suc_list, [])})
        final_df = final_df.append(succseed_df)
    return final_df


def crop_rot_figures():
    if not os.path.exists('./figures/'):
        # Create the temp directory if it does not exist
        os.makedirs('./figures/')
    # Assuming in_path_1 and in_path_2 are your file paths
    inital = pd.read_csv('./' + out_path + '/crop_rot_freq_init.csv')
    optimized = pd.read_csv('./' + out_path + '/crop_rot_freq_' + str(tolerance) + '.csv')

    inital['opt'] = "initial"
    optimized['opt'] = "optimized"

    # Concatenate the dataframes
    test = pd.concat([inital, optimized], ignore_index=True)

    # Get unique crop types and pre_suc values
    crop_types = test['crop_t_from'].unique()
    pre_suc_values = test['pre_suc'].unique()
    # Create subplots
    for croptype in crop_types:
        fig, axes = plt.subplots(1, len(pre_suc_values), figsize=(15, 10))
        plt.subplots_adjust(bottom=0.2)
        for j, pre_suc_value in enumerate(pre_suc_values):
            ax = axes[j]
            #ax.figure()
            ax.hist([test[(test['crop_t_from'] == croptype) & (test['opt'] == 'initial') & (test['pre_suc'] == pre_suc_value)]['value'],
                      test[(test['crop_t_from'] == croptype) & (test['opt'] == 'optimized') & (test['pre_suc'] == pre_suc_value)]['value']],
                     bins=np.arange(len(test['value'].unique()) + 1),
                     alpha=0.7, label=['initial', 'optimized'], align='left')

            #ax.set_xticks(np.arange(len(test['value'].unique())))
            ax.set_xlabel('', fontsize=14)  # Adjust x-axis label font size
            ax.set_ylabel('Frequency', fontsize=14)  # Adjust y-axis label font size
            ax.set_xticklabels(test[(test['crop_t_from'] == croptype) &
                                    #(test['opt'] == 'initial') &
                                    (test['pre_suc'] == pre_suc_value)]['value'].unique(), rotation=45, ha="right", fontsize=14)
            #ax.set_xlabel('Value')
            ax.set_title(f"{str(croptype)} - {pre_suc_value}", fontsize=16)
            ax.legend()
        #plt.show()
        plt.savefig('./figures/' + str(croptype) + '.png')
        # Adjust layout to prevent overlapping
        #plt.tight_layout()
        #plt.show()
