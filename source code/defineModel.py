import tensorflow as tf
tf.compat.v1.enable_eager_execution()
from Hyperpar import *
from scipy.signal import gaussian
import numpy as np
"""
All tensor objects will be defined as global variables.
Any source can access to these values for restoring or
saving or training tensors.
"""
def tensorSummaries(tensor,name):
    with tf.name_scope(name):
        mean = tf.reduce_mean(tensor)
        tf.summary.scalar('mean', mean)
        with tf.name_scope('stddev'):
            stddev = tf.sqrt(tf.reduce_mean(tf.square(tensor - mean)))
        tf.summary.scalar('stddev', stddev)
        tf.summary.scalar('max', tf.reduce_max(tensor))
        tf.summary.scalar('min', tf.reduce_min(tensor))
        tf.summary.histogram('histogram', tensor)

def wrapTensorToTensorboard():
    pass

############# Loss & Dropout probability  ##############################
is_train_step = tf.Variable(initial_value=False, trainable=False, dtype=tf.bool)

loss_weight = tf.Variable(initial_value=False, trainable=False, dtype=tf.float32)
class_threshold = tf.Variable(initial_value=False, trainable=False, dtype=tf.float32)
p_dropout = tf.Variable(initial_value=False, trainable=False, dtype=tf.float32)

############################### INPUT ##################################
input_data_train = tf.Variable(initial_value=tf.zeros((batch_size, target_size, 1), dtype=tf.float32), trainable=False, name="trainingData")
input_data_eval = tf.Variable(initial_value=tf.zeros((batch_size, target_size, 1), dtype=tf.float32), trainable=False, name="testData")

label_data_train = tf.Variable(initial_value=tf.zeros((batch_size, 1, target_size), dtype=tf.float32), trainable=False, name="labelData")
label_data_eval = tf.Variable(initial_value=tf.zeros((batch_size, 1, target_size), dtype=tf.float32), trainable=False, name="labelDataEval")

input_ref_data_train = tf.Variable(initial_value=tf.zeros((batch_size, target_size, 1), dtype=tf.float32), trainable=False, name="TrainRefData")
input_ref_data_eval = tf.Variable(initial_value=tf.zeros((batch_size, target_size, 1), dtype=tf.float32), trainable=False, name="TestRefData")

smoothing_filter = tf.constant(gaussian(filter_size_b,50)/np.sum(gaussian(filter_size_b,50)),tf.float32, shape=[filter_size_b,1,1], name='smoothing_filter')

###################### STEM FOR REFGENEDEPTH ###########################
conv1_ref_weight = tf.Variable(tf.compat.v1.get_variable("Conv_REF_1", shape=[3, 1, conv1_ref_features], initializer= tf.keras.initializers.VarianceScaling(scale=1.0, mode='fan_avg',distribution='uniform')))


############################ STEM Layer  ###############################
conv1_weight = tf.Variable(tf.compat.v1.get_variable("Conv_STEM1", shape=[5, 1, conv1_features], initializer= tf.keras.initializers.VarianceScaling(scale=1.0, mode='fan_avg',distribution='uniform')))
conv2_weight = tf.Variable(tf.compat.v1.get_variable("Conv_STEM2", shape=[3, conv1_features, conv2_features], initializer= tf.keras.initializers.VarianceScaling(scale=1.0, mode='fan_avg',distribution='uniform')))
conv1_bias = tf.Variable(tf.compat.v1.get_variable("Conv_STEM1_bias", shape=[conv1_features], initializer= tf.keras.initializers.VarianceScaling(scale=1.0, mode='fan_avg',distribution='uniform')))
conv2_bias = tf.Variable(tf.compat.v1.get_variable("Conv_STEM2_bias", shape=[conv2_features], initializer= tf.keras.initializers.VarianceScaling(scale=1.0, mode='fan_avg',distribution='uniform')))

############################ Inception 1 ###############################
layer1_width = conv2_features + conv1_ref_features
conv1a_weight = tf.Variable(tf.compat.v1.get_variable("Conv_1A", shape=[4, layer1_width, conv1a_features], initializer= tf.keras.initializers.VarianceScaling(scale=1.0, mode='fan_avg',distribution='uniform')))
conv1b_weight = tf.Variable(tf.compat.v1.get_variable("Conv_1B", shape=[3, layer1_width, conv1b_features], initializer= tf.keras.initializers.VarianceScaling(scale=1.0, mode='fan_avg',distribution='uniform')))
conv1c_weight = tf.Variable(tf.compat.v1.get_variable("Conv_1C", shape=[2, layer1_width, conv1c_features], initializer= tf.keras.initializers.VarianceScaling(scale=1.0, mode='fan_avg',distribution='uniform')))
convMax1_weight = tf.Variable(tf.compat.v1.get_variable("Conv_max_W1", shape=[2, layer1_width, convMax1_features], initializer= tf.keras.initializers.VarianceScaling(scale=1.0, mode='fan_avg',distribution='uniform')))
convAvg1_weight = tf.Variable(tf.compat.v1.get_variable("Conv_avg_W1", shape=[2, layer1_width, convAvg1_features], initializer= tf.keras.initializers.VarianceScaling(scale=1.0, mode='fan_avg',distribution='uniform')))

conv1a_bias = tf.Variable(tf.compat.v1.get_variable("conv1a_bias", shape=[conv1a_features], initializer= tf.keras.initializers.VarianceScaling(scale=1.0, mode='fan_avg',distribution='uniform')))
conv1b_bias = tf.Variable(tf.compat.v1.get_variable("conv2b_bias", shape=[conv1b_features], initializer= tf.keras.initializers.VarianceScaling(scale=1.0, mode='fan_avg',distribution='uniform')))
conv1c_bias = tf.Variable(tf.compat.v1.get_variable("conv3c_bias", shape=[conv1c_features], initializer= tf.keras.initializers.VarianceScaling(scale=1.0, mode='fan_avg',distribution='uniform')))

############################ Inception 2 ###############################
layer2_width = conv1a_features + conv1b_features + conv1c_features + (layer1_width*2)
conv2a_weight = tf.Variable(tf.compat.v1.get_variable("Conv_2A", shape=[4, layer2_width, conv2a_features], initializer= tf.keras.initializers.VarianceScaling(scale=1.0, mode='fan_avg',distribution='uniform')))
conv2b_weight = tf.Variable(tf.compat.v1.get_variable("Conv_2B", shape=[3, layer2_width, conv2b_features], initializer= tf.keras.initializers.VarianceScaling(scale=1.0, mode='fan_avg',distribution='uniform')))
conv2c_weight = tf.Variable(tf.compat.v1.get_variable("Conv_2C", shape=[2, layer2_width, conv2c_features], initializer= tf.keras.initializers.VarianceScaling(scale=1.0, mode='fan_avg',distribution='uniform')))
convMax2_weight = tf.Variable(tf.compat.v1.get_variable("Conv_max_W2", shape=[2, layer2_width, convMax2_features], initializer= tf.keras.initializers.VarianceScaling(scale=1.0, mode='fan_avg',distribution='uniform')))
convAvg2_weight = tf.Variable(tf.compat.v1.get_variable("Conv_avg_W2", shape=[2, layer2_width, convAvg2_features], initializer= tf.keras.initializers.VarianceScaling(scale=1.0, mode='fan_avg',distribution='uniform')))

############################ Inception 3 ###############################
layer3_width = conv2a_features + conv2b_features + (layer2_width *2 ) + conv2c_features
conv3a_weight = tf.Variable(tf.compat.v1.get_variable("Conv_3A", shape=[5, layer3_width, conv3a_features],initializer= tf.keras.initializers.VarianceScaling(scale=1.0, mode='fan_avg',distribution='uniform')))
conv3b_weight = tf.Variable(tf.compat.v1.get_variable("Conv_3B", shape=[4, layer3_width, conv3b_features],initializer= tf.keras.initializers.VarianceScaling(scale=1.0, mode='fan_avg',distribution='uniform')))
conv3c_weight = tf.Variable(tf.compat.v1.get_variable("Conv_3C", shape=[3, layer3_width, conv3c_features],initializer= tf.keras.initializers.VarianceScaling(scale=1.0, mode='fan_avg',distribution='uniform')))

############################ Inception 4 ###############################
layer4_width = conv3a_features + conv3b_features + conv3c_features + layer3_width
conv4a_weight = tf.Variable(tf.compat.v1.get_variable("Conv_4A", shape=[4, layer4_width, conv4a_features], initializer= tf.keras.initializers.VarianceScaling(scale=1.0, mode='fan_avg',distribution='uniform')))
conv4b_weight = tf.Variable(tf.compat.v1.get_variable("Conv_4B", shape=[3, layer4_width, conv4b_features], initializer= tf.keras.initializers.VarianceScaling(scale=1.0, mode='fan_avg',distribution='uniform')))
conv4c_weight = tf.Variable(tf.compat.v1.get_variable("Conv_4C", shape=[2, layer4_width, conv4c_features],initializer= tf.keras.initializers.VarianceScaling(scale=1.0, mode='fan_avg',distribution='uniform')))

############################ Inception 5 ###############################
layer5_width = conv4a_features + conv4b_features + conv4c_features +layer4_width
conv5a_weight = tf.Variable(tf.compat.v1.get_variable("Conv_5A", shape=[4, layer5_width, conv5a_features],initializer= tf.keras.initializers.VarianceScaling(scale=1.0, mode='fan_avg',distribution='uniform')))
conv5b_weight = tf.Variable(tf.compat.v1.get_variable("Conv_5B", shape=[3, layer5_width, conv5b_features],initializer= tf.keras.initializers.VarianceScaling(scale=1.0, mode='fan_avg',distribution='uniform')))
conv5c_weight = tf.Variable(tf.compat.v1.get_variable("Conv_5C", shape=[2, layer5_width, conv5c_features],initializer= tf.keras.initializers.VarianceScaling(scale=1.0, mode='fan_avg',distribution='uniform')))

############################ Inception 6  ###############################
layer6_width = conv5a_features + conv5b_features + conv5c_features +layer5_width
conv6a_weight = tf.Variable(tf.compat.v1.get_variable("Conv_6A", shape=[10, layer6_width, conv6a_features],initializer= tf.keras.initializers.VarianceScaling(scale=1.0, mode='fan_avg',distribution='uniform')))
conv6b_weight = tf.Variable(tf.compat.v1.get_variable("Conv_6B", shape=[8, layer6_width, conv6b_features],initializer= tf.keras.initializers.VarianceScaling(scale=1.0, mode='fan_avg',distribution='uniform')))
conv6c_weight = tf.Variable(tf.compat.v1.get_variable("Conv_6C", shape=[6, layer6_width, conv6c_features],initializer= tf.keras.initializers.VarianceScaling(scale=1.0, mode='fan_avg',distribution='uniform')))
conv6d_weight = tf.Variable(tf.compat.v1.get_variable("Conv_6D", shape=[5, layer6_width, conv6d_features],initializer= tf.keras.initializers.VarianceScaling(scale=1.0, mode='fan_avg',distribution='uniform')))

############################ Inception 6  ###############################
layer7_width = conv6a_features + conv6b_features + conv6c_features + conv6d_features
conv7a_weight = tf.Variable(tf.compat.v1.get_variable("Conv_7A", shape=[10, layer7_width, conv7a_features],initializer= tf.keras.initializers.VarianceScaling(scale=1.0, mode='fan_avg',distribution='uniform')))
conv7b_weight = tf.Variable(tf.compat.v1.get_variable("Conv_7B", shape=[8, layer7_width, conv7b_features],initializer= tf.keras.initializers.VarianceScaling(scale=1.0, mode='fan_avg',distribution='uniform')))
conv7c_weight = tf.Variable(tf.compat.v1.get_variable("Conv_7C", shape=[6, layer7_width, conv7c_features],initializer= tf.keras.initializers.VarianceScaling(scale=1.0, mode='fan_avg',distribution='uniform')))
conv7d_weight = tf.Variable(tf.compat.v1.get_variable("Conv_7D", shape=[5, layer7_width, conv7d_features],initializer= tf.keras.initializers.VarianceScaling(scale=1.0, mode='fan_avg',distribution='uniform')))


############################ Fully Connected ###############################
layer_full_width = conv7a_features + conv7b_features + conv7c_features + conv7d_features#conv6a_features + conv6b_features + conv6c_features + conv6d_features
resulting_width = 1

final_conv_size = resulting_width * ( layer_full_width )
full1_weight = tf.Variable(tf.compat.v1.get_variable("Full_W1", shape=[final_conv_size, fully_connected_size1],initializer= tf.keras.initializers.VarianceScaling(scale=1.0, mode='fan_avg',distribution='uniform')))
full2_weight = tf.Variable(tf.compat.v1.get_variable("Full_W2", shape=[fully_connected_size1, fully_connected_size2], initializer= tf.keras.initializers.VarianceScaling(scale=1.0, mode='fan_avg',distribution='uniform')))

############################ Output ###############################
output_weight = tf.Variable(tf.compat.v1.get_variable("Full_Output", shape=[fully_connected_size2, threshold_division], initializer= tf.keras.initializers.VarianceScaling(scale=1.0, mode='fan_avg',distribution='uniform')))
output_bias = tf.Variable(tf.random.truncated_normal([threshold_division], stddev=0.1, dtype=tf.float32))



def peakPredictConvModel(input_data_depth, input_data_ref, test=False, smoothing=True, normalize=False, eps=0.00001):
    """
    Define structure of convolution model.
    It will return two final output tensor.
     1. Sub output tensor for training steps
     2. Final output tensor for training and test steps
    :param logger:
    :param input_data:
    :return: Tensor of the output layer
    """
    input_data_depth = input_data_depth + init_depth 

    if smoothing:
        depth_tensor_max = tf.nn.pool(input_data_depth, [filter_size_a], strides=[1], padding='SAME', pooling_type='MAX')
        depth_tensor_smooth = tf.nn.conv1d(depth_tensor_max, smoothing_filter, stride=1, padding='SAME')
        input_data_depth = depth_tensor_smooth

    if normalize:
        depth_mean, depth_var = tf.nn.moments(input_data_depth, [1])
        input_data_depth = (input_data_depth -  depth_mean)/tf.sqrt(depth_var + eps)
    

    #Stem of read depth data
    conv1 = tf.nn.conv1d(input_data_depth, conv1_weight, stride=1, padding='SAME')
    conv1 = tf.keras.layers.BatchNormalization(axis=-1, momentum=0.9, epsilon=1e-5, center=True, scale=True,
                                       trainable=is_train_step)(conv1)
    conv1 = tf.nn.relu(conv1)
    print("Stem 1 : {}".format(conv1.shape))

    conv2 = tf.nn.conv1d(conv1, conv2_weight, stride=1, padding='SAME')
    conv2 = tf.keras.layers.BatchNormalization(axis=-1, momentum=0.9, epsilon=1e-5, center=True, scale=True,
                                       trainable=is_train_step)(conv2)
    conv2 = tf.nn.relu(conv2)
    max_pool1 = tf.nn.pool(conv2, [max_pool_size_stem], strides=[max_pool_size_stem],
            padding='SAME', pooling_type='MAX')
    print("Stem 2 : {}".format(max_pool1.shape))

    #Stem of ref gene data
    conv1_ref = tf.nn.conv1d(input_data_ref, conv1_ref_weight, stride=1, padding='SAME')
    conv1_ref = tf.keras.layers.BatchNormalization(axis=-1, momentum=0.9, epsilon=1e-5, center=True, scale=True,
                                       trainable=is_train_step)(conv1_ref)
    conv1_ref = tf.nn.relu(conv1_ref)
    max_pool1_ref = tf.nn.pool(conv1_ref, [max_pool_size_stem], strides=[max_pool_size_stem],
            padding='SAME', pooling_type='MAX')
    print("Stem_Ref 1 : {}".format(max_pool1_ref.shape))

    input_concat = tf.concat([max_pool1, max_pool1_ref],axis = 2) #max_pool1
    print("Stem_concat : {}".format(input_concat.shape))

    # Inception modules 1 to 6
    concat1 = concatLayer_A(input_concat, conv1a_weight, convMax1_weight, conv1b_weight, convAvg1_weight, conv1c_weight, 2)

    concat2 = concatLayer_A(concat1, conv2a_weight, convMax2_weight, conv2b_weight, convAvg2_weight, conv2c_weight, 2)

    concat3 = concatLayer_B(concat2, conv3a_weight, conv3b_weight, conv3c_weight, 3)
    concat3 = concat3 + tf.nn.pool(tf.concat([concat1 for x in range(4)], axis=2), [6], strides=[6], pooling_type='AVG', padding='SAME')

    concat4 = concatLayer_B(concat3, conv4a_weight, conv4b_weight, conv4c_weight, 2)

    concat5 = concatLayer_B(concat4, conv5a_weight, conv5b_weight, conv5c_weight, 2)
    concat5 = concat5 + tf.nn.pool(tf.concat([concat3 for x in range(3)], axis=2), [4], strides=[4], pooling_type='AVG', padding='SAME')

    concat6 = concatLayer_C(concat5, conv6a_weight, conv6b_weight, conv6c_weight, conv6d_weight, 5)

    concat7 = concatLayer_C(concat6, conv7a_weight, conv7b_weight, conv7c_weight, conv7d_weight, 5)
    concat7 = concat7 + tf.nn.pool(tf.concat([concat5 for x in range(2)], axis=2), [25], strides=[25], pooling_type='AVG', padding='SAME')

    concat7 = tf.nn.pool(concat7, [5], strides=[5], padding='SAME', pooling_type='AVG')
    print("Final Avg pool : {}".format(concat7.shape))

    final_conv_shape = concat7.get_shape().as_list()
    final_shape = final_conv_shape[1] * final_conv_shape[2]
    flat_output = tf.reshape(concat7, [final_conv_shape[0] , final_shape])

    fully_connected1 = tf.matmul(flat_output, full1_weight)
    fully_connected1 = tf.keras.layers.BatchNormalization(axis=-1, momentum=0.9, epsilon=1e-5, center=True, scale=True,
                                       trainable=is_train_step)(fully_connected1)
    fully_connected1 = tf.nn.relu(fully_connected1)
    print("Fully connected A :{}".format(fully_connected1.shape))

    fully_connected2 = tf.matmul(fully_connected1, full2_weight)
    fully_connected2 = tf.keras.layers.BatchNormalization(axis=-1, momentum=0.9, epsilon=1e-5, center=True, scale=True,
                                       trainable=is_train_step)(fully_connected2)
    fully_connected2 = tf.nn.relu(fully_connected2)
    print("Fully connected B :{}".format(fully_connected2.shape))

    final_threshold_output = tf.matmul(fully_connected2, output_weight)
    final_threshold_output = tf.add(final_threshold_output, output_bias)
    print("Threshold Output :{}".format(final_threshold_output.shape))

    return final_threshold_output


def concatLayer_B(source_layer, conv1_w, conv2_w, conv3_w, pooling_size):
    """
    Define concat layer which like Inception module.
    :param source_layer:
    :param conv1_w:
    :param conv2_w:
    :param conv1_b:
    :param conv2_b:
    :param pooling_size:
    :return:
    """
    conv1 = tf.nn.conv1d(source_layer, conv1_w, stride=pooling_size, padding='SAME')
    conv1 = tf.keras.layers.BatchNormalization(axis=-1, momentum=0.9, epsilon=1e-5, center=True, scale=True,
                                       trainable=is_train_step)(conv1)
    conv1 = tf.nn.relu(conv1)

    conv2 = tf.nn.conv1d(source_layer, conv2_w, stride=pooling_size, padding='SAME')
    conv2 = tf.keras.layers.BatchNormalization(axis=-1, momentum=0.9, epsilon=1e-5, center=True, scale=True,
                                       trainable=is_train_step)(conv2)
    conv2 = tf.nn.relu(conv2)

    conv3 = tf.nn.conv1d(source_layer, conv3_w, stride=pooling_size, padding='SAME')
    conv3 = tf.keras.layers.BatchNormalization(axis=-1, momentum=0.9, epsilon=1e-5, center=True, scale=True,
                                       trainable=is_train_step)(conv3)
    conv3 = tf.nn.relu(conv3)

    max_pool = tf.nn.pool(source_layer, [pooling_size], strides=[pooling_size],
            padding='SAME', pooling_type='MAX')

    concat = tf.concat([conv1, conv2, conv3,  max_pool], axis=2)
    print("Concat Type A :{}".format(concat.shape))
    return concat


def concatLayer_C(source_layer, conv1_w, conv2_w, conv3_w, conv4_w, pooling_size):
    """
    Define concat layer which like Inception module.
    :param source_layer:
    :param conv1_w:
    :param conv_max_w:
    :param conv2_w:
    :param conv_avg_w:
    :param conv1_b:
    :param conv_max_b:
    :param conv2_b:
    :param conv_avg_b:
    :param pooling_size:
    :return:
    """
    conv1 = tf.nn.conv1d(source_layer, conv1_w, stride=pooling_size, padding='SAME')
    conv1 = tf.keras.layers.BatchNormalization(axis=-1, momentum=0.9, epsilon=1e-5, center=True, scale=True,
                                       trainable=is_train_step)(conv1)
    conv1 = tf.nn.relu(conv1)

    conv2 = tf.nn.conv1d(source_layer, conv2_w, stride=pooling_size, padding='SAME')
    conv2 = tf.keras.layers.BatchNormalization(axis=-1, momentum=0.9, epsilon=1e-5, center=True, scale=True,
                                       trainable=is_train_step)(conv2)
    conv2 = tf.nn.relu(conv2)

    conv3 = tf.nn.conv1d(source_layer, conv3_w, stride=pooling_size, padding='SAME')
    conv3 = tf.keras.layers.BatchNormalization(axis=-1, momentum=0.9, epsilon=1e-5, center=True, scale=True,
                                       trainable=is_train_step)(conv3)
    conv3 = tf.nn.relu(conv3)

    conv4 = tf.nn.conv1d(source_layer, conv4_w, stride=pooling_size, padding='SAME')
    conv4 = tf.keras.layers.BatchNormalization(axis=-1, momentum=0.9, epsilon=1e-5, center=True, scale=True,
                                       trainable=is_train_step)(conv4)
    conv4 = tf.nn.relu(conv4)

    concat = tf.concat([conv1, conv2, conv3, conv4], axis=2)
    print("Concat Type C :{}".format(concat.shape))
    return concat


def concatLayer_A(source_layer, conv1_w, conv_max_w, conv2_w, conv_avg_w, conv3_w, pooling_size):
    """
    Define concat layer which like Inception module.
    :param source_layer:
    :param pooling_size:
    :return:
    """
    conv1 = tf.nn.conv1d(source_layer, conv1_w, stride=2, padding='SAME')
    conv1 = tf.keras.layers.BatchNormalization(axis=-1, momentum=0.9, epsilon=1e-5, center=True, scale=True,
                                       trainable=is_train_step)(conv1)
    conv1 = tf.nn.relu(conv1)

    conv2 = tf.nn.conv1d(source_layer, conv2_w, stride=2, padding='SAME')
    conv2 = tf.keras.layers.BatchNormalization(axis=-1, momentum=0.9, epsilon=1e-5, center=True, scale=True,
                                       trainable=is_train_step)(conv2)
    conv2 = tf.nn.relu(conv2)

    conv3 = tf.nn.conv1d(source_layer, conv3_w, stride=2, padding='SAME')
    conv3 = tf.keras.layers.BatchNormalization(axis=-1, momentum=0.9, epsilon=1e-5, center=True, scale=True,
                                       trainable=is_train_step)(conv3)
    conv3 = tf.nn.relu(conv3)

    max_pool = tf.nn.pool(source_layer, [pooling_size], strides=[pooling_size,],
            padding='SAME', pooling_type='MAX')

    avg_pool = tf.nn.pool(source_layer, [pooling_size], strides=[pooling_size,],
            padding='SAME', pooling_type='AVG')

    concat = tf.concat([conv1, avg_pool, conv2, max_pool, conv3], axis=2)
    print("Concat Type A :{}".format(concat.shape))
    return concat


def generateOutput(threshold_tensor, depth_tensor, div=10, input_size=12000, batch_size_in=batch_size, smoothing=True, normalize=False, t=1, eps=0.00001):
    """
    It generate
    :param threshold_tensor: This tensor represents read-depth thresholds which have size of 'div'
    :param depth_tensor: This tensor represents
    :param div:
    :return:
    """
    depth_tensor = depth_tensor + init_depth 
    
    if smoothing:
        depth_tensor_max = tf.nn.pool(depth_tensor, [filter_size_a], strides=[1], padding='SAME', pooling_type='MAX')
        depth_tensor_smooth = tf.nn.conv1d(depth_tensor_max, smoothing_filter, stride=1, padding='SAME')
        depth_tensor = depth_tensor_smooth

    if normalize:
        depth_mean, depth_var = tf.nn.moments(depth_tensor, [1])
        depth_tensor = (depth_tensor -  depth_mean)/tf.sqrt(depth_var + eps)
        #depth_max = tf.reduce_max(depth_tensor)
        #depth_tensor = ((depth_tensor - depth_mean)/(depth_max + 0.00001))*100
        
    depth_tensor = tf.reshape(depth_tensor,[batch_size_in ,div, input_size//div])
    threshold_tensor = tf.reshape(threshold_tensor,[batch_size_in,div,1])
    print("depth tensor :{}".format(depth_tensor.shape))
    print("threshold tensor :{}".format(threshold_tensor.shape))

    result_tensor = tf.subtract(depth_tensor, threshold_tensor, name="results")
    print("result tensor before:{}".format(result_tensor.shape))
    result_tensor = tf.reshape(result_tensor,[batch_size_in, 1, input_size])
    print("result tensor after:{}\n".format(result_tensor.shape))

    if not (t == 1):
        result_tensor = result_tensor * t 

    return result_tensor



def aggregatedLoss(y_true, y_pred):
    # Calculate binary cross-entropy
    loss_a = tf.reduce_mean(tf.keras.losses.binary_crossentropy(y_true, y_pred))

    return loss_a


######################## Tensor graph for training steps #################################
# ...

# Tensor graph for training steps
model_output = peakPredictConvModel(input_data_train, input_ref_data_train)
prediction_before_sigmoid = generateOutput(model_output, input_data_train, div=threshold_division, smoothing=True)
prediction = tf.nn.sigmoid(prediction_before_sigmoid)
loss = aggregatedLoss(label_data_train, prediction_before_sigmoid)

# Get trainable variables
trainable_variables = tf.compat.v1.trainable_variables()

# Get update_ops for batch normalization
update_ops = tf.compat.v1.get_collection(tf.compat.v1.GraphKeys.UPDATE_OPS)

# Ensure update_ops are run before applying gradients
with tf.control_dependencies(update_ops):
    # Use a tape to compute gradients
    with tf.GradientTape() as tape:
        loss = aggregatedLoss(label_data_train, prediction_before_sigmoid)
        # loss = topKAverageLoss(label_data_train, prediction_before_sigmoid)
    # Get gradients
    gradients = tape.gradient(loss, trainable_variables)

    # Define optimizer and apply gradients
    optimizer = tf.keras.optimizers.Adam(learning_rate)
    train_step = optimizer.apply_gradients(zip(gradients, trainable_variables))

####################### Tensor graph for test steps ######################################
test_model_output = peakPredictConvModel(input_data_eval, input_ref_data_eval)
test_prediction = tf.nn.sigmoid(generateOutput(test_model_output, input_data_eval, div=threshold_division, smoothing=True))




