import matplotlib.pyplot as plt
import numpy as np

from __functions import *


def no_x_after_y_check(vals_, crop_x, crop_y):
    values_x = np.where(np.array(vals_) == crop_x, 1, 0)
    values_y = np.where(np.array(vals_) == crop_y, 1, 0)

    n_years = len(vals_)
    for year_ in range(n_years):
        if year_ <= n_years - 2:
            if values_y[year_] > 0.5 and values_x[year_ + 1] > 0.5:
                return 1
    return 0


def min_return_t_for_x_check(vals_, t, x):
    # we do t-1 because then setting t is then more intuitive. Setting t = 2 means that a crop can be grown every
    # second year [1, 0, 1]
    values_x = np.where(np.array(vals_) == x, 1, 0)

    n_years = len(vals_)
    t = t - 1
    # adapt t to n_years, so we don't index more than the list length
    t_adapted_to_n_years = t
    for year_ in range(n_years):
        while t_adapted_to_n_years + year_ > n_years-1:
            t_adapted_to_n_years = t_adapted_to_n_years - 1

        if values_x[year_] > 0.5 and any(values_x[year_ + add] > 0.5 for add in range(1, t_adapted_to_n_years + 1)):
            # calculates the shortest brake between the same crop; i.e. the minimum return time
            return shortest_sequence(values_x.astype(bool))[0]+1

    return t+1


def longest_seq_x(vals_, x):
    crop_ = x
    values = np.where(np.array(vals_) == crop_, 1, 0)
    longest_seq, max_start_index = longest_sequence(values)
    return longest_seq


def check_CropRotRules(historic_croptypes_dict):
    # here we check for each field if our predefined crop rotation rules were violated;
    # i need a function that returns a dict in the following structure:
    # constr_is_enforcable_dict {field_1: [1, 1, 1], field_2: [1, 0, 1],...}
    # for each field we check if each of our e.g. three constraints is violated by the farmer; if yes -> if no -> 0;
    # this dict will be used with the lazyconstraints; if dict 1 -> don't enforce this constraint
    print('checking rules')
    print(historic_croptypes_dict)
    rapeseed_minret_violation_list = []
    potato_minret_violation_list = []
    beet_minret_violation_list = []
    legume_minret_violation_list = []
    sunflower_minret_violation_list = []
    sunflower_rp_violation_list = []
    no_rapeseed_after_maize_list = []
    no_beets_after_rapeseed = []
    no_legume_after_rapeseed_list = []
    no_sunflowers_after_legumes_list = []
    no_rapeseed_after_potatoes_list = []
    no_maize_after_beets_list = []




    cereals_longest_seq_list = []
    maize_longest_seq_list = []
    springCereal_longest_seq_list = []
    for key, value in historic_croptypes_dict.items():
        rapeseed_minret_violation_list.append(min_return_t_for_x_check(value, t=3, x=crop_names_dict['rapeseed']))

        potato_minret_violation_list.append(min_return_t_for_x_check(value, t=4, x=crop_names_dict['potato']))

        beet_minret_violation_list.append(min_return_t_for_x_check(value, t=4, x=crop_names_dict['beets']))

        legume_minret_violation_list.append(min_return_t_for_x_check(value, t=2, x=crop_names_dict['legumes']))

        sunflower_minret_violation_list.append(min_return_t_for_x_check(value, t=4, x=crop_names_dict['sunflowers']))

        # no rapeseed after sunflowers and vice versa
        if no_x_after_y_check(value, crop_x=crop_names_dict['sunflowers'], crop_y=crop_names_dict['rapeseed']) == 1 or \
                no_x_after_y_check(value, crop_x=crop_names_dict['rapeseed'], crop_y=crop_names_dict['sunflowers']) == 1:
            sunflower_rp_violation_list.append(1)
        else:
            sunflower_rp_violation_list.append(0)

        # no rapeseed after potato and vice versa
        if no_x_after_y_check(value, crop_x=crop_names_dict['potato'], crop_y=crop_names_dict['rapeseed']) == 1 or \
                no_x_after_y_check(value, crop_x=crop_names_dict['rapeseed'], crop_y=crop_names_dict['potato']) == 1:
            no_rapeseed_after_potatoes_list.append(1)
        else:
            no_rapeseed_after_potatoes_list.append(0)

        # no rapeseed after sunflowers and vice versa
        if no_x_after_y_check(value, crop_x=crop_names_dict['beets'], crop_y=crop_names_dict['rapeseed']) == 1 or \
                no_x_after_y_check(value, crop_x=crop_names_dict['rapeseed'], crop_y=crop_names_dict['beets']) == 1:
            no_beets_after_rapeseed.append(1)
        else:
            no_beets_after_rapeseed.append(0)

        # no maize after beets and vice versa
        if no_x_after_y_check(value, crop_x=crop_names_dict['beets'], crop_y=crop_names_dict['maize']) == 1 or \
                no_x_after_y_check(value, crop_x=crop_names_dict['maize'], crop_y=crop_names_dict['beets']) == 1:
                no_maize_after_beets_list.append(1)
        else:
            no_maize_after_beets_list.append(0)

        no_rapeseed_after_maize_list.append(no_x_after_y_check(value, crop_x=crop_names_dict['rapeseed'], crop_y=crop_names_dict['maize']))
        no_legume_after_rapeseed_list.append(no_x_after_y_check(value, crop_x=crop_names_dict['legumes'], crop_y=crop_names_dict['rapeseed']))
        no_sunflowers_after_legumes_list.append(no_x_after_y_check(value, crop_x=crop_names_dict['sunflowers'], crop_y=crop_names_dict['legumes']))
        cereals_longest_seq_list.append(longest_seq_x(value, x=crop_names_dict['winter_cereals']))
        maize_longest_seq_list.append(longest_seq_x(value, x=crop_names_dict['maize']))
        springCereal_longest_seq_list.append(longest_seq_x(value, x=crop_names_dict['spring_cereals']))

    print(sum([1 if value != 3 else 0 for value in rapeseed_minret_violation_list]), 'violations of the rapeseed constraint')
    print(sum([1 if value != 4 else 0 for value in potato_minret_violation_list]), 'violations of the potato constraint')
    print(sum([1 if value != 4 else 0 for value in beet_minret_violation_list]), 'violations of the beet constraint')
    print(sum([1 if value != 2 else 0 for value in legume_minret_violation_list]), 'violations of the legumes constraint')
    print(sum([1 if value != 4 else 0 for value in sunflower_minret_violation_list]), 'violations of the sunflower min return constraint')
    print(sum(sunflower_rp_violation_list), 'violations of the sunflower-rapeseed constraint')
    print(sum(no_beets_after_rapeseed), 'violations of the no_beets_after_rapeseed')
    print(sum(no_rapeseed_after_maize_list), 'no_rapeseed_after_maize_list')
    print(sum(no_legume_after_rapeseed_list), 'no_legume_after_rapeseed_list')
    print(sum(no_sunflowers_after_legumes_list), 'no_sunflowers_after_legumes_list')
    print(sum(no_rapeseed_after_potatoes_list), 'no_rapeseed_after_potatoes_list')

    violation_dict = {id_: (val1, val2, val3, val4, val5, val6, val7, val8, val9, val10, val11, val12) for id_, val1, val2, val3, val4, val5, val6, val7, val8, val9, val10, val11, val12 in zip(list(historic_croptypes_dict.keys()), rapeseed_minret_violation_list, potato_minret_violation_list,
                                                                beet_minret_violation_list, legume_minret_violation_list, sunflower_rp_violation_list, sunflower_minret_violation_list,
                                                                                                                   no_rapeseed_after_maize_list, no_legume_after_rapeseed_list, no_sunflowers_after_legumes_list,
                                                                                                                                                                     no_beets_after_rapeseed, no_rapeseed_after_potatoes_list, no_maize_after_beets_list)}
    longest_seq_dict = {id_: (val1, val2, val3) for id_, val1, val2, val3 in zip(list(historic_croptypes_dict.keys()),
                                                                     maize_longest_seq_list, cereals_longest_seq_list, springCereal_longest_seq_list)}
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
