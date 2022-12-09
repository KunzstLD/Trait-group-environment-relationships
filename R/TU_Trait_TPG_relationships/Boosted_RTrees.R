# https://xgboost.readthedocs.io/en/latest/tutorials/model.html
# https://github.com/mlr-org/mlr3
# XGBOOST, boosted regression trees
mlr_learners$get("regr.xgboost")
learner <- lrn("regr.xgboost")
