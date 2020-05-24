
How to use the Model object
===========================

This page describes the usage of the ``Model`` object. This object provides a
pure Python API, so some basic knowledge about Python is assumed.

Model paradigm
--------------

A ``Model`` consists one or more *features*. For instance, an *aircraft* model
might have the feature *propulsion* and the feature *wing*. Such a feature
consists of one or more *properties*. Concrete values can be assigned to these
properties. So, in our example propulsion might have the properties *thrust*
and *type*, and the *wing* feature might have properties *span* and
*mount_point* (e.g. for an engine, landing gear, or some other technical
system). We may assign some values to each of these properties. For instance,
we may set thrust to :math:`50\,\textrm{kN}` or the span to
:math:`20\,\textrm{m}`. How do we set these values in a Python script?

.. code:: python

    # Get an instance of the Model object
    model = Model()

    # From the model instance, we can create a 'propulsion' instance. We can
    # now set actual value to the properties 'trust' and propulsion 'type'.
    prop = model.set_feature('propulsion')
    prop.set('thrust', 50e3)
    prop.set('type', 'turbofan engine')

    # Now we create a new 'wing'. We can set the 'span' and add mount points
    wing = model.add_feature('wing')
    wing.set('span', 20)
    wing.add('mount_point', (4, 2, 4))
    wing.add('mount_point', (4, 6, 4))
    wing.add('mount_point', (4, 20, 4))


.. figure:: https://raw.githubusercontent.com/airinnova/model-framework/master/src/mframework/ressources/aircraft_model.svg
   :width: 500
   :alt: Simple aircraft model
   :align: center

   Simple aircraft model with a *propulsion* features, and a *wing* feature

Notice the method ``add_feature()``. When calling this method, we create a
*instance* of the feature we want to add (here ``'propulsion'`` and
``'wing'``). The feature instance has ``set()`` and ``add()`` methods to assign
a value to some property. The model object will check the value. That means
that you must provide a number for the property *trust*. If you try to assign,
say a string (``'50e3'``), an error will be thrown, and the model will not
continue to be built.


.. figure:: https://raw.githubusercontent.com/airinnova/model-framework/master/src/mframework/ressources/model_api_hierarchy.svg
   :alt: Model hierarchy
   :align: center

   A model object can have multiple features, and each feature can have
   multiple properties

Model and feature methods
-------------------------

The ``Model`` object and its features provide only few methods. The main thing
to remember is that there are ``set*`` and ``add*`` methods. The ``set*``
method will always apply if there can only be one instance, and the ``add*``
method will apply if there can be multiple instances. For example, an aircraft
can have more than one wing, hence the ``add_feature()`` method applies. To
create another wing, we simply call the method again. However, some feature may
only exists once. In our dummy aircraft model, we impose the restriction of
only adding *one* propulsion feature (arguably, an aircraft model could have
multiple propulsion instances, but here we assume otherwise). To highlight that
propulsion is a *singleton* feature, we use the ``set_feature()`` method.
Trying to use ``add_feature('propulsion')`` will result in an error and the
model will not continue to build. Property values in a feature can be assigned
with the ``add()`` and ``set()`` methods, depending on whether the properties
are singleton or not. In the example above, the wing may have multiple mount
points, but there can only be a single wing span.

+---------------+----------------------------------+---------------------------------+
|               | **Model**                        | **Feature**                     |
+---------------+----------------------------------+---------------------------------+
| Singleton     | ``set_feature('feature_name')``  | ``set('property_name', value)`` |
+---------------+----------------------------------+---------------------------------+
| Non-singleton | ``add_feature('feature_name')``  | ``add('property_name', value)`` |
+---------------+----------------------------------+---------------------------------+

You can retrieve feature instances from your model at any time. Note that you
must have created these instances first with the methods discussed above. The
table below shows the available methods you may use.

..
    TODO
    * More detailed explanation of difference between 'get()' and 'iter()'

+---------------+--------------------------+---------------------------+
|               | **Model**                | **Feature**               |
+---------------+--------------------------+---------------------------+
| Singleton     | ``get('feature_name')``  | ``get('property_name')``  |
+---------------+--------------------------+---------------------------+
| Non-singleton | ``iter('feature_name')`` | ``iter('property_name')`` |
+---------------+--------------------------+---------------------------+

Running the model
-----------------

The last thing you need to know is how to run the model. Once you set up the
entire, you can call the ``run()`` method. This method will start the actual
evaluation of the model.

.. code:: python

    model = Model()
    ...  # Add features and assign values here
    results = model.run()

Results
-------

The ``run()`` method returns an object with which you can interact pretty much
in the same way as the ``Model`` object. Results are group into "features"
which have properties. You can retrieve data using the ``get()`` and ``iter()``
methods mentioned above.



