import numpy as np
# from sklearn.metrics import auc, accuracy_score, precision_score, recall_score, precision_recall_curve, auc

def evaluate_performance(dec, labels, type='classification'):
    if type == 'classification':
        precision = precision_score(labels, np.where(dec >= 0, 1, -1))
        recall = recall_score(labels, np.where(dec >= 0, 1, -1))
        sensitivity = recall
        specificity = specificity_score(labels, np.where(dec >= 0, 1, -1))
        fscore = 2 * (precision * recall) / (precision + recall)
        bac = (sensitivity + specificity) / 2
        # auc_score = auc(labels, dec)
        # accuracy = accuracy_score(labels, np.where(dec >= 0, 1, -1))
        aupr = aupr_score(labels, dec)
        # ci = concordance_index(dec, labels)
        auc_score = 0
        accuracy = 0
    else:
        precision, recall, sensitivity, specificity, fscore, bac, auc_score, accuracy, aupr = 0, 0, 0, 0, 0, 0, 0, 0, 0
        # ci = concordance_index(dec, labels)

    return {'precision': precision, 'recall': recall, 'sensitivity': sensitivity, 'specificity': specificity, 'fscore': fscore, 'bac': bac, 'auc': auc_score, 'accuracy': accuracy, 'aupr': aupr, }


def precision_score(y_true, y_pred):
    tp = np.sum((y_true == 1) & (y_pred >= 0))
    tp_fp = np.sum(y_pred >= 0)
    if tp_fp == 0:
        print('warning: No positive predict label.')
        ret = 0
    else:
        ret = tp / tp_fp
    return ret


def recall_score(y_true, y_pred):
    tp = np.sum((y_true == 1) & (y_pred >= 0))
    tp_fn = np.sum(y_true == 1)
    if tp_fn == 0:
        print('warning: No postive true label.')
        ret = 0
    else:
        ret = tp / tp_fn
    return ret


def specificity_score(y_true, y_pred):
    tn = np.sum((y_true == -1) & (y_pred < 0))
    tn_fp = np.sum(y_true == -1)
    if tn_fp == 0:
        print('warning: No negative true label.')
        specificity = 0
    else:
        specificity = tn / tn_fp
    return specificity


def fscore_score(y_true, y_pred):
    tp = np.sum((y_true == 1) & (y_pred >= 0))
    tp_fp = np.sum(y_pred >= 0)
    tp_fn = np.sum(y_true == 1)
    if tp_fp == 0:
        print('warning: No positive predict label.')
        precision = 0
    else:
        precision = tp / tp_fp
    if tp_fn == 0:
        print('warning: No postive true label.')
        recall = 0
    else:
        recall = tp / tp_fn
    if precision + recall == 0:
        print('warning: precision + recall = 0.')
        ret = 0
    else:
        ret = 2 * precision * recall / (precision + recall)
    return ret


# def aupr(dec, label):
#     precision, recall, _ = precision_recall_curve(label, dec)
#     return auc(recall, precision)


def aupr_score(y_true, y_score):
    # sort in descending order of predicted scores
    desc_score_indices = np.argsort(y_score)[::-1]
    y_score = y_score[desc_score_indices]
    y_true = y_true[desc_score_indices]

    # calculate precision and recall for each threshold
    precisions = []
    recalls = []
    num_positives = np.sum(y_true)
    true_positives = 0
    print(type(y_true),type(y_score),type(y_true[0]))
    for i in range(len(y_score)):
        if y_true[i] == 1:
            true_positives += 1
        precisions.append(true_positives / (i + 1))
        recalls.append(true_positives / num_positives)

    # calculate AUPR using trapezoidal rule
    aupr = 0
    for i in range(1, len(recalls)):
        aupr += (recalls[i] - recalls[i-1]) * ((precisions[i] + precisions[i-1]) / 2)

    return aupr
