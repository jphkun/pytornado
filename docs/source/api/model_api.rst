Model
=====

Below you will find a comprehensive list of all
available features and properties. The model object has the following features:



.. mermaid::

    graph TD
    A[Model]
    A --> F0[wing]


Feature: wing
-------------

.. image:: https://raw.githubusercontent.com/airinnova/model-framework/master/src/mframework/ressources/icons/notes.svg
   :align: left
   :alt: description

*Description*: Add a wing

.. image:: https://raw.githubusercontent.com/airinnova/model-framework/master/src/mframework/ressources/icons/point.svg
   :align: left
   :alt: singleton

*Singleton*: False

.. image:: https://raw.githubusercontent.com/airinnova/model-framework/master/src/mframework/ressources/icons/lifebuoy.svg
   :align: left
   :alt: required

*Required*: True

Property: segment_vertices
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. mermaid::

    graph LR
    A[Model]
    A --> F1[wing] 
    F1 --> P1[segment_vertices] 


.. image:: https://raw.githubusercontent.com/airinnova/model-framework/master/src/mframework/ressources/icons/notes.svg
   :align: left
   :alt: description

*Description*: Add a wing segment

.. image:: https://raw.githubusercontent.com/airinnova/model-framework/master/src/mframework/ressources/icons/point.svg
   :align: left
   :alt: singleton

*Singleton*: False

.. image:: https://raw.githubusercontent.com/airinnova/model-framework/master/src/mframework/ressources/icons/lifebuoy.svg
   :align: left
   :alt: required

*Required*: True

.. image:: https://raw.githubusercontent.com/airinnova/model-framework/master/src/mframework/ressources/icons/clipboard-check.svg
   :align: left
   :alt: schema

*Schema*:

======== ==============
**type** <class 'dict'>
======== ==============

