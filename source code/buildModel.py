import os
#os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
import shutil
import pandas as pd
import numpy as np
import tensorflow as tf
tf.compat.v1.enable_eager_execution()
import matplotlib.pyplot as plt
import random
from scipy.ndimage import maximum_filter1d
from scipy.signal import gaussian
from utilities import extractChrClass  
from defineModel import *
from Hyperpar import *

def run(dir_name, logger, num_grid=0, K_fold_in=10, cross_valid = False, label_fix_num=None):
    """
    This is a main function to build convolution neural network model
    for peak prediction.
    :param dir_name:
    :param logger:
    :param num_grid:
    :param K_fold_in:
    :return:
    """
    PATH = os.path.abspath(dir_name)

    dir_list = os.listdir(PATH)

    for dir in dir_list:
        True

    input_list = {}
    for dir in dir_list:
        dir = PATH + '/' + dir
        input_list[dir] = extractChrClass(dir)

    print("smoothing filter size A :{}".format(filter_size_a))
    print("smoothing filter size B :{}".format(filter_size_b))
    print("Learning Rate :{}".format(learning_rate))
    print("Threshold division :{}".format(threshold_division))

    ###################### Define Session and Initialize Graph ###############################
    init_glob = tf.global_variables_initializer()
    init_local = tf.local_variables_initializer()
    sess = tf.Session()
    sess.run(init_glob)
    sess.run(init_local)

    ###################### Data Load and Split for Cross validation ###### ###################
    input_data_names = []
    ref_data_names = []
    label_data_names = []

    for dir in input_list:
        label_num = 0
        for chr in input_list[dir]:
            for cls in input_list[dir][chr]:
                label_num += 1
                input_file_name = "{}/{}_{}_grid{}.ct".format(dir, chr, cls, num_grid)
                ref_file_name = "{}/ref_{}_{}_grid{}.ref".format(dir, chr, cls, num_grid)
                label_file_name = "{}/label_{}_{}_grid{}.lb".format(dir, chr, cls, num_grid)
                input_data_names.append(pd.read_csv(input_file_name))
                ref_data_names.append(pd.read_csv(ref_file_name))
                label_data_names.append(pd.read_csv(label_file_name))
        logger.info("DIRECTORY (TARGET) [{}]# of labels : <{}>".format(label_num,dir))

    if label_fix_num != None:
        step = 0
        fin = len(input_data_names)
        while fin-step > label_fix_num:
            rand_idx = random.randint(0, fin - step - 1)
            input_data_names.pop(rand_idx)
            ref_data_names.pop(rand_idx)
            label_data_names.pop(rand_idx)
            step += 1
            
        
    K_fold = K_fold_in
    input_data_list, label_data_list, ref_data_list = splitTrainingData(input_data_names, label_data_names
            , ref_data_names, Kfold=K_fold)
    if not os.path.isdir(os.getcwd() + "/models"):      ### Make directory for save the model
        os.mkdir(os.getcwd() + "/models")

    #################### Training start with cross validation ################################
    for i in range(K_fold):

        training_data = []
        training_ref = []
        training_label = []
        test_data = []
        test_ref = []
        test_label = []
        for j in range(K_fold):
            if i == j:
                test_data += input_data_list[j]
                test_ref += ref_data_list[j]
                test_label += label_data_list[j]
            else:
                training_data += input_data_list[j]
                training_ref += ref_data_list[j]
                training_label += label_data_list[j]

        if not os.path.isdir(os.getcwd() + "/models/model_{}".format(i)):
            os.mkdir(os.getcwd() + "/models/model_{}".format(i))

        logger.info("[{}]# of labels for training, [{}]# of labels for testing.".format(len(training_data), len(test_data)))
        training(sess, loss, prediction, test_prediction, train_step, training_data, training_label
                , training_ref, test_data, test_label, test_ref, logger, num_grid, i)

        if cross_valid == False:
            break


def training(sess, loss, prediction, test_prediction, train_step, train_data_list, train_label_list
        , train_ref_list, test_data_list, test_label_list, test_ref_list, logger, num_grid, step_num):
    """
    :param train_data_list:
    :param train_label_list:
    :param test_data_list:
    :param test_label_list:
    :param train_step:
    :param loss:
    :param prediction:
    :param test_prediction:
    :param logger:
    :param num_grid:
    :param step_num:
    :return:
    """

    train_loss = []
    train_spec = []
    train_sens = []
    test_spec = []
    test_sens = []

    loss_containor_for_mean = []
    spec_containor_for_mean = []
    sens_containor_for_mean = []

    test_spec_containor_for_mean = []
    test_sens_containor_for_mean = []

    saver = tf.train.Saver(max_to_keep=10)
    
    #pnRate_sum = 0
    #for i in range(len(train_data_list)):
    #    rand_y = []
    #    rand_y.append(np.repeat(train_label_list[i][['peak']].values.transpose(),5))
    #    rand_y = np.array(rand_y).reshape(label_data_train.shape)
    #    pnRate_sum += pnRate(rand_y)
    #pnRate_mean = pnRate_sum/len(train_data_list)
    #logger.info("average PN rate : {}".format(pnRate_mean))
    
    # Start of the training process
    for i in range(generations):
        rand_index = np.random.choice(len(train_data_list), size=batch_size)

        rand_x = []
        rand_ref = []
        rand_y = []

        for j in range(batch_size):
            rand_x.append(train_data_list[rand_index[j]]['readCount'].values)
            rand_ref.append(train_ref_list[rand_index[j]]['refGeneCount'].values)
            rand_y.append(np.repeat(train_label_list[rand_index[j]][['peak']].values.transpose(),5))

        rand_x = np.array(rand_x).reshape(input_data_train.shape)
        rand_ref = np.array(rand_ref).reshape(input_ref_data_train.shape)
        rand_y = np.array(rand_y).reshape(label_data_train.shape)

        train_dict = {input_data_train: rand_x, label_data_train: rand_y, input_ref_data_train: rand_ref,
                loss_weight:pnRate(rand_y), is_train_step:True, p_dropout:0.6}
                #loss_weight:pnRate_mean, is_train_step:True, p_dropout:0.6}

        sess.run([train_step], feed_dict=train_dict)

        temp_train_loss, temp_train_preds = sess.run([loss, prediction], feed_dict=train_dict)
        temp_train_stat = getStat(temp_train_preds, rand_y, num_grid=num_grid)

        loss_containor_for_mean.append(temp_train_loss)
        spec_containor_for_mean.append(temp_train_stat['spec'])
        if temp_train_stat['sens'] != -1:
            sens_containor_for_mean.append(temp_train_stat['sens'])

        # Recording results of test data
        eval_index = np.random.choice(len(test_data_list), size=batch_size)

        eval_x = []
        eval_ref = []
        eval_y = []

        for j in range(batch_size):
            eval_x.append(test_data_list[eval_index[j]]['readCount'].values)
            eval_ref.append(test_ref_list[eval_index[j]]['refGeneCount'].values)
            eval_y.append(np.repeat(test_label_list[eval_index[j]][['peak']].values.transpose(),5))

        eval_x = np.array(eval_x).reshape(input_data_eval.shape)
        eval_ref = np.array(eval_ref).reshape(input_ref_data_eval.shape)
        eval_y = np.array(eval_y).reshape(label_data_eval.shape)

        test_dict = {input_data_eval: eval_x, label_data_eval: eval_y, input_ref_data_eval: eval_ref
                , is_train_step:False}

        test_preds = sess.run(test_prediction, feed_dict=test_dict)
        temp_test_stat = getStat(test_preds, eval_y, num_grid=num_grid)
        
        
        test_spec_containor_for_mean.append(temp_test_stat['spec'])
        if temp_test_stat['sens'] != -1:
            test_sens_containor_for_mean.append(temp_test_stat['sens'])

        if (i + 1) % eval_every == 0:
            if len(sens_containor_for_mean) == 0:
                sens_mean = -1.
            else:
                sens_mean = sum(sens_containor_for_mean)/float(len(sens_containor_for_mean))
                loss_mean = sum(loss_containor_for_mean)/float(len(loss_containor_for_mean))
            spec_mean = sum(spec_containor_for_mean)/float(len(spec_containor_for_mean))

            if len(test_sens_containor_for_mean) == 0:
                test_sens_mean = -1.
            else:
                test_sens_mean = sum(test_sens_containor_for_mean)/float(len(test_sens_containor_for_mean))
            test_spec_mean = sum(test_spec_containor_for_mean)/float(len(test_spec_containor_for_mean))

            logger.info('Generation # {}. Loss: {:.2f}. Test: SENS:{:.2f} SPEC:{:.2f}| Train: SENS:{:.2f} SPEC:{:.2f}\n'.format(i+1, loss_mean, test_sens_mean, test_spec_mean, sens_mean, spec_mean))

            train_loss.append(loss_mean)
            train_spec.append(spec_mean)
            train_sens.append(sens_mean)
            test_spec.append(test_spec_mean)
            test_sens.append(test_sens_mean)

            loss_containor_for_mean.clear()
            spec_containor_for_mean.clear()
            sens_containor_for_mean.clear()
            test_spec_containor_for_mean.clear()
            test_sens_containor_for_mean.clear()
        if i % 100 == 0 and i != 1:
            save_path = saver.save(sess, os.getcwd() + "/models/model{}.ckpt".format(step_num,step_num), global_step=i)

    visualizeTrainingProcess(eval_every, generations, test_sens, test_spec, train_sens, train_spec
            , train_loss, K_fold=str(step_num))
    visualizePeakResult(batch_size, input_data_eval, num_grid, label_data_eval, sess, test_data_list
            , test_label_list, test_ref_list, test_prediction, k=len(test_data_list), K_fold=str(step_num))

    logger.info("Saving CNNPeaks")
    save_path = saver.save(sess, os.getcwd() + "/models/model{}.ckpt".format(step_num,step_num))
    logger.info("Model saved in path : %s" % save_path)


def getTensorStat(logits, targets, batch_size_in=batch_size):
    """
    Return accuracy of the result.
    Acc = ( TP + TN ) / ( TP + TN + FN + FP )
    ( TP + TN + FN + FP ) = num_grid
    :param logits:
    :param targets:
    :return:
    """

    threshold_tensor = tf.constant(class_threshold)

    sensDummy = tf.constant(-1.,dtype=tf.float64)

    for i in range(batch_size_in):
        TP = tf.count_nonzero(tf.cast(tf.logical_and(tf.less_equal(threshold_tensor, logits[i][0]),
                tf.less_equal(threshold_tensor, targets[i][0])), dtype=tf.float32))
        TN = tf.count_nonzero(tf.cast(tf.logical_and(tf.less_equal(logits[i][0],threshold_tensor),
                tf.less_equal(targets[i][0], threshold_tensor)), dtype=tf.float32))
        FN = tf.count_nonzero(tf.cast(tf.logical_and(tf.less(logits[i][0], threshold_tensor),
                tf.less(threshold_tensor, targets[i][0])), dtype=tf.float32))
        FP = tf.count_nonzero(tf.cast(tf.logical_and(tf.less(threshold_tensor, logits[i][0]),
                tf.less(targets[i][0], threshold_tensor)), dtype=tf.float32))

    sens = tf.cond( tf.equal(tf.add(TP,FN),0) , true_fn=lambda:sensDummy, false_fn=lambda:tf.divide(TP,tf.add(TP,FN)))
    spec = tf.divide(TN,tf.add(TN,FP))

    return sens, spec


def getStat(logits, targets, batch_size_in=batch_size, num_grid=0):
    """
    Return accuracy of the result.
    Acc = ( TP + TN ) / ( TP + TN + FN + FP )
    ( TP + TN + FN + FP ) = num_grid
    :param logits:
    :param targets:
    :return:
    """
    logits = logits.reshape(batch_size_in, num_grid)
    targets = targets.reshape(batch_size_in, num_grid)

    TP = 0.
    TN = 0.
    FN = 0.
    FP = 0.

    for i in range(batch_size_in):
        for index in range(len(logits[0])):
            if (logits[i][index]) >= class_threshold and targets[i][index] >= class_threshold:
                TP += 1
            elif (logits[i][index]) >= class_threshold and targets[i][index] < class_threshold:
                FP += 1
            elif (logits[i][index]) < class_threshold and targets[i][index] >= class_threshold:
                FN += 1
            elif (logits[i][index]) < class_threshold and targets[i][index] < class_threshold:
                TN += 1
            else:
                pass

    if TP+FN == 0:
        return {'sens': -1, 'spec': TN/(TN+FP), 'acc':(TP+TN)/(TP+TN+FN+FP)}
    else:
        return {'sens': TP/(TP+FN), 'spec': TN/(TN+FP), 'acc':(TP+TN)/(TP+TN+FN+FP)}


def pnRate(targets, batch_size_in=batch_size):
    """
    Return the The ratio of Negative#/ Positive#.
    It will be used for weights of loss function to adjust
    between sensitivity and specificity.
    :param targets:
    :param num_grid:
    :return:
    """
    positive_count = 0.

    for i in range(batch_size_in):
        for index in range(len(targets[i][0])):
            if targets[i][0][index] > 0.5:
                positive_count += 1

    # For the label only has negative samples.
    if positive_count == 0.:
        return 1

    negative_count = len(targets[0][0])*batch_size_in - positive_count

    ### TODO :: adjust these SHITTY EQUATION.
    return (negative_count / positive_count)


def classValueFilter(output_value):
    """
    For output of model, probabilities of a final vector will be changed
    as binary values by checking whether elements of vector are higher or lower than
    class_threshold that defined in hyperparameters.py.
    :param output_value:
    :return: a binary vector that indicates having peak or not.
    """

    class_value_list = []

    for index in range(output_value.shape[2]):
        if output_value[0][0][index] >= class_threshold:
            class_value_list.append(1)
        elif output_value[0][0][index] < class_threshold:
            class_value_list.append(0)

    return class_value_list


def splitTrainingData(data_list, label_list, ref_list, Kfold=15):
    """
    If Kfold is zero, it just split two parts
    :param list_data:
    :param Kfold:
    :return:
    """
    print("################## THE NUMBER OF LABEL DATA : {}".format(len(data_list)))

    size = len(data_list)
    counter = size / Kfold

    test_data = []
    test_label = []
    test_ref = []

    for i in range(Kfold - 1):
        test_data_temp = []
        test_ref_temp = []
        test_label_temp = []
        while True:
            if counter <= 0:
                test_ref.append(test_ref_temp)
                test_data.append(test_data_temp)
                test_label.append(test_label_temp)
                counter = size // Kfold
                break

            pop_index = random.randint(0, len(data_list)-1)
            test_ref_temp.append(ref_list.pop(pop_index))
            test_data_temp.append(data_list.pop(pop_index))
            test_label_temp.append(label_list.pop(pop_index))
            counter -= 1

    test_data.append(data_list)
    test_ref.append(ref_list)
    test_label.append(label_list)

    return test_data, test_label, test_ref


def patternVis(reads, labels, refs, preds, num_grid=12000,min_peak_size=10, max_peak_num=50):
    reads = reads.reshape(num_grid)
    refs = refs.reshape(num_grid)
    reads = maximum_filter1d(reads, filter_size_a)  ## MAX POOL to extract boarder lines
    smoothing_filter = gaussian(filter_size_b, 50) / np.sum(gaussian(filter_size_b, 50))
    reads = np.convolve(reads, smoothing_filter, mode='same')  ## Smoothing boarder lines
    
    #plt.subplot(2,1,2, autoscale_on=True,position=[0,2,1,0.2])
    #plt.plot(preds.reshape(num_grid))
    
    plt.subplot(2,1,2, autoscale_on=True,position=[0.05,2,0.9,0.2])
    plt.imshow(preds.reshape(num_grid)[np.newaxis,:], cmap="jet", vmax=1, vmin=0, aspect="auto")
            
    preds = classValueFilter(preds)
    labels = classValueFilter(labels)

    ############# Peak post processing ##########
    peak_num = 0

    if preds[0] > 0:
        peak_size = 1
        peak_num += 1
    else:
        peak_size = 0

    for pred_index in range(len(preds)-1):
        if preds[pred_index+1] > 0:
            peak_size += 1
        else:
            if peak_size < min_peak_size:
                for j in range(peak_size):
                    preds[pred_index-j] = 0
                peak_size=0
            else:
                peak_num += 1
                peak_size = 0

    if peak_num > max_peak_num:
        for j in range(len(preds)):
            preds[j] = 0
    #############################################

    y_index = []
    y = []
    for index in range(len(preds)):
        if preds[index] > 0:
            y_index.append(index)
            y.append(labels[index])

    ref_index = []
    ref = []
    for index in range(len(refs)):
        if refs[index] > 0:
            ref_index.append(index)
            ref.append(-1)
            
    plt.subplot(2,1,1,position=[0,2.2,1,1])
    plt.plot(reads.reshape(num_grid).tolist(),'k', markersize=2, linewidth=1)
    plt.plot(y_index,y, 'r.', label='Model prediction')
    plt.plot(ref_index, ref, 'b|', markersize=8)

    onPositive = False
    start = 0
    end = 0
    for j in range(len(labels)):
        if labels[j] == 1 and not onPositive:
            start = j
            onPositive = True
        elif (labels[j] == 0 or j == len(labels)-1) and onPositive:
            end = j
            onPositive = False
            plt.axvspan(start, end, color='red', alpha=0.3)


    plt.title('Peak prediction result by regions')
    plt.xlabel('Regions')
    plt.ylabel('Read Count')
    plt.legend(loc='upper right')
    plt.show()
    plt.clf()


def visualizePeakResult(batch_size, input_data_eval, num_grid, label_data_eval, sess, test_data_list, test_label_list,
                        test_ref_list, test_prediction, k = 1, K_fold="", min_peak_size=10, max_peak_num=50):
    """
    :param batch_size:
    :param input_data_eval:
    :param num_grid:
    :param label_data_eval:
    :param sess:
    :param test_data_list:
    :param test_label_list:
    :param test_prediction:
    :param k:
    :return:
    """

    if k > 0:
        for i in range(k):
            show_x = test_data_list[i]['readCount'].values
            show_x = show_x.reshape(input_data_eval.shape)

            show_ref = test_ref_list[i]['refGeneCount'].values
            show_ref = show_ref.reshape(input_data_eval.shape)

            show_y = test_label_list[i][['peak']].values.transpose()
            show_y = np.repeat(show_y, 5)
            show_y = show_y.reshape(label_data_train.shape)
            show_dict = {input_data_eval: show_x, input_ref_data_eval: show_ref, label_data_eval: show_y,
                         is_train_step: False, p_dropout:1}
            show_preds = sess.run(test_prediction, feed_dict=show_dict)

            show_x = show_x.reshape(num_grid)
            smoothing_filter = gaussian(filter_size_b, 50) / np.sum(gaussian(filter_size_b, 50))
            show_x = maximum_filter1d(show_x, filter_size_a)  ## MAX POOL to extract boarder lines
            show_x = np.convolve(show_x, smoothing_filter, mode='same')  ## Smoothing boarder lines

            show_ref = show_ref.reshape(num_grid)
            
            plt.subplot(2,1,2, autoscale_on=True,position=[0,2,1,0.2])
            plt.imshow(show_preds.reshape(num_grid)[np.newaxis,:], cmap="jet", vmax=1, vmin=0, aspect="auto")
            
            show_preds = classValueFilter(show_preds)
            show_y = classValueFilter(show_y)

            ############# Peak post processing ##########
            peak_num = 0

            if show_preds[0] > 0:
                peak_size = 1
                peak_num += 1
            else:
                peak_size = 0

            for pred_index in range(len(show_preds)-1):
                if show_preds[pred_index+1] > 0:
                    peak_size += 1
                else:
                    if peak_size < min_peak_size:
                        for j in range(peak_size):
                            show_preds[pred_index-j] = 0
                        peak_size=0
                    else:
                        peak_num += 1
                        peak_size = 0

            if peak_num > max_peak_num:
                for j in range(len(show_preds)):
                    show_preds[j] = 0
            #############################################

            y_index = []
            y = []
            for index in range(len(show_preds)):
                if show_preds[index] > 0:
                    y_index.append(index)
                    y.append(show_y[index])

            ref_index = []
            ref = []
            for index in range(len(show_ref)):
                if show_ref[index] > 0:
                    ref_index.append(index)
                    ref.append(-1)
            
            plt.subplot(2,1,1,position=[0,2.2,1,1])
            plt.plot(show_x.reshape(num_grid).tolist(),'k', markersize=2, linewidth=1)
            plt.plot(y_index,y, 'r.', label='Model prediction')
            plt.plot(ref_index, ref, 'b|', markersize=8)

            onPositive = False
            start = 0
            end = 0
            for j in range(len(show_y)):
                if show_y[j] == 1 and not onPositive:
                    start = j
                    onPositive = True
                elif (show_y[j] == 0 or j == len(show_y)-1) and onPositive:
                    end = j
                    onPositive = False
                    plt.axvspan(start, end, color='red', alpha=0.3)


            plt.title('Peak prediction result by regions')
            plt.xlabel('Regions')
            plt.ylabel('Read Count')
            plt.legend(loc='upper right')
            plt.savefig('models/model_{}/peak{}.png'.format(K_fold,i), dpi=400)
            plt.show()
            plt.clf()


def visualizeTrainingProcess(eval_every, generations, test_sens, test_spec, train_sens, train_spec, train_loss, K_fold =""):
    """
    Create matplotlib figures about a plot of loss function values and accuracy values.
    :param eval_every:
    :param generations:
    :param test_acc:
    :param train_acc:
    :param train_loss:
    :return:
    """
    eval_indices = range(0, generations, eval_every)

    for i in range(len(train_sens) - 1):
        if train_sens[i+1] == -1:
            train_sens[i+1] = train_sens[i]
        if test_sens[i+1] == -1:
            test_sens[i+1] = test_sens[i]

    plt.plot(eval_indices, train_loss, 'k-')
    plt.title('Cross entropy Loss per generation')
    plt.xlabel('Generation')
    plt.ylabel('Loss')
    plt.savefig('models/model_{}/LossPerGen.png'.format(K_fold), dpi=400)
    plt.show()
    plt.clf()

    plt.plot(eval_indices, test_sens, label='Test Set sensitivity')
    plt.plot(eval_indices, train_sens, label='Train Set sensitivity')
    plt.title('Train and Test Sensitivity')
    plt.xlabel('Generation')
    plt.ylabel('Sensitivity')
    plt.legend(loc='lower right')
    plt.ylim([0,1])
    plt.savefig('models/model_{}/SensPerGen.png'.format(K_fold), dpi=400)
    plt.show()
    plt.clf()

    plt.plot(eval_indices, test_spec, label='Test Set specificity')
    plt.plot(eval_indices, train_spec, label='Train Set specificity')
    plt.title('Train and Test specificity')
    plt.xlabel('Generation')
    plt.ylabel('Specificity')
    plt.legend(loc='lower right')
    plt.ylim([0,1])
    plt.savefig('models/model_{}/SpecPerGen.png'.format(K_fold), dpi=400)
    plt.show()
    plt.clf()

    plt.plot(eval_indices, [2*(test_sens[i]*test_spec[i])/(test_sens[i]+test_spec[i]) for i in range(len(test_sens))], label='Test Set F1 Score')
    plt.plot(eval_indices, [2*(train_sens[i]*train_spec[i])/(train_sens[i]+train_spec[i]) for i in range(len(train_sens))], label='Train Set F1 Score')
    plt.title('F1 Score')
    plt.xlabel('Generation')
    plt.ylabel('F1 Score')
    plt.legend(loc='lower right')
    plt.ylim([0,1])
    plt.savefig('models/model_{}/F1scorePerGen.png'.format(K_fold), dpi=400)
    plt.show()
    plt.clf()
